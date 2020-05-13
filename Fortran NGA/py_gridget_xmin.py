"""
source: Fortran code "gridget_1min.f" found at:
https://earth-info.nga.mil/GandG///wgs84/gravitymod/egm2008/egm08_wgs84.html

I translated it into python, and then adapted the code to what I wanted it
to do. 

c-----------------------------------------------------------------------
c     ORIGINAL PROGRAM:                           SIMON HOLMES, JUL 2007
c     MODIFIED FOR CORNER-CELL REGISTRATION       SIMON HOLMES, MAY 2008
c     MUTLI-OUTPUT OPTION                         SIMON HOLMES, MAY 2008
#     TRANSLATION INTO PYTHON              XAVIER DE LABRIOLLE, MAY 2020
c-----------------------------------------------------------------------
"""
import numpy as np
from scipy.io import FortranFile


# =============================================================================
# MAIN
# =============================================================================
line_5000 = "#" +"-" * 60
print("\n")
print(line_5000)
print("\t Extract a grid of pre-computed geoid undulations")
print(line_5000)

# FILE PARAMETERS
path_in  = "../../Fortran NGA/Fortran grid"
n_in     = "Und_min1x1_egm2008_isw=82_WGS84_TideFree_SE"
path_out = "../Rendered/temp"
nout     = "pyOUTPUT.txt"

# UNFORMATTED GRID PARAMETERS
nrows = 10801 
ncols = 21600
dlat = 1/60 # degrees
dlon = 1/60 # degrees

# PROMPT FOR GEOGRAPHIC REGION
print("Enter limits (degrees): west east south north")
ok = False
while not(ok):
    print("longitude - latitude: -180 +180 -90 90")
    dwest, deast, dsouth, dnorth = map( float, input().split(' ') )
    
    if ( (dnorth>90) or (dsouth<-90) or 
        (dwest>deast) or (dsouth>dnorth) or 
        (dwest<-180) or (deast>180)):
        print("\t INVALID, cartopy cannot plot beyond 180th meridian")
    else:
        ok = True
print("\n")

# PROMPT FOR RESOLUTION
print ("Enter sample spacing : dlat_out dlon_out")
ok = False
while not(ok):
    print("Spacing in minutes : 1 - 600(10 degrees)")
    dlat_out, dlon_out = map( float, input().split(' ') )
    
    if ((dlat_out<=0) or (dlat_out>600) or
        (dlon_out<=0) or (dlon_out>600) or
        (int(dlat_out)!=dlat_out) or
        (int(dlon_out)!=dlon_out) ):
        print("\t INVALID, enter integers within (1, 600)")
    else: 
        ok = True
print("\n")

# convert to minutes
dlat_out = dlat_out / 60
dlon_out = dlon_out / 60

# =============================================================================
# INITIALISE DIMENSIONS OF EXTRACTED GRID
# =============================================================================
north_i = round((90-dnorth) / dlat)
south_i = round((90-dsouth) / dlat) # included
west_j  = round((180+dwest) / dlon)
east_j  = round((180+deast) / dlon) # included

if (west_j < 0):     west_j += ncols
if (west_j > ncols): west_j -= ncols

if (east_j < 0):     east_j += ncols
if (east_j > ncols): east_j -= ncols

north_m = 90 - north_i*dlat
south_m = 90 - south_i*dlat
west_m  = west_j*dlon - 180
east_m  = east_j*dlon - 180

print(west_j, east_j, south_i, north_i)
print(west_m, east_m, south_m, north_m)

flat = np.arange(north_m, south_m-dlat/3, -dlat_out/(dlat*60) )
flon = np.arange(west_m,  east_m+dlat/3,   dlon_out/(dlon*60) )
if (west_j > east_j):
    flon_a = np.arange(west_m, 180, dlon_out/(dlon*60))
    flon_b = np.arange(-180, east_m+dlat/3,dlon_out/(dlon*60))
    flon = [*flon_a,*flon_b]

irow  = len(flat)
jcol  = len(flon)
shape = (irow, jcol)

fetch = lambda jj : west_j + jj*round(dlat_out/dlat) + int(ncols/2)
skip_lon = round(dlon_out/dlon) - 1

print(line_5000)
print("Geometry of extracted grid:\n")
print(  f"Latitude of northern boundary = {flat[0]:.3f} (Degrees)",
      f"\nLatitude of southern boundary = {flat[-1]:.3f} (Degrees)",
      f"\nLongitude of western boundary = {flon[0]:.3f} (Degrees)",
      f"\nLongitude of eastern boundary = {flon[-1]:.3f} (Degrees)",
      f"\n\n           Latitude spacing = {dlat_out*60} (Minutes)",
        f"\n          Longitude spacing = {dlon_out*60} (Minutes)",
      f"\n\nGrid is {irow} rows x {jcol} columns of values")


# =============================================================================
# OPEN INPUT AND OUTPUT FILES
# =============================================================================
fnul1  = f"{path_in}/{n_in}"
file_1 = FortranFile(fnul1, 'r')

fnul10 = f"{path_out}/{nout}"
file_10 = open(fnul10, "w+") # output file

#print(line_5000)
#print(f"Input Sequential Binary Data File : \n\t{n_in}\n")
#print(f"Output Extracted Ascii Data File  : \n\t{nout}\n")

90
# =============================================================================
# EXTRACT GRID ONE ROW AT A TIME
# =============================================================================
print(line_5000)
print("Extracting ...\n")
print("Grid corner values:")
print("LAT\t\tLONG\t\tUND")

# LOOP-READ THROUGH UNINTERESTING ROWS DATA 
for i in range(0, north_i):
    _ = file_1.read_record(dtype=np.float32)


# LOOP-EXTRACT INTERESTING ROW DATA
for ii in range (0, irow):  
    # This line reads one whole line of the unformatted binary file
    grid = file_1.read_record(dtype=np.float32)
     
    temp = np.zeros(jcol)
    for jj in range (0, jcol):
        fetch_j=fetch(jj)
        temp[jj] = grid[fetch_j%ncols]
    
    for j in range (0, jcol):
        file_10.write(f"{flat[ii]:.6f}\t{flon[j]:.6f}\t{temp[j]:.6f}\n")
        if ( ((j==0) or (j==jcol-1)) and ((ii==0) or (ii==irow-1)) ): 
            print(f"{flat[ii]:.6f}\t{flon[j]:.6f}\t{temp[j]:.6f}")
    
    # LOOP READ UNWANTED ROWS
    if (ii != irow-1):
        for s in range (0, skip_lon):
            _ = file_1.read_record(dtype=np.float32)

file_1.close()
file_10.close()

"""
print("\n")
print(line_5000)
print("This is where I would be showing stats but I didn't code it.",
      "\nAlso I dont see the point.\n")

print(line_5000)
print("\t Normal Termination")
print(line_5000)
print("Credtis: Xavier de Labriolle")
"""
