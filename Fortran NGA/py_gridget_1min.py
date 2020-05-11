"""
source: Fortran code for gridget 1min

c-----------------------------------------------------------------------
c     ORIGINAL PROGRAM:                           SIMON HOLMES, JUL 2007
c     MODIFIED FOR CORNER-CELL REGISTRATION       SIMON HOLMES, MAY 2008
c     MUTLI-OUTPUT OPTION                         SIMON HOLMES, MAY 2008
#     TRANSLATION INTO PYTHON              XAVIER DE LABRIOLLE, MAY 2020
c-----------------------------------------------------------------------

"""
import numpy as np
from scipy.io import FortranFile

term_width = 50
line_5000 = "-" * term_width

def stats():
    pass

# =============================================================================
# MAIN
# =============================================================================
print(line_5000)
print("\t Execution")
print(line_5000)

# FILE PARAMETERS
path_in  = './Fortran grid'
n_in     = 'Und_min1x1_egm2008_isw=82_WGS84_TideFree_SE'
path_out = '.'
nout     = 'pyOUTPUT.DAT'

# UNFORMATTED GRID PARAMETERS
nrows = 10801 
ncols = 21600
dlat = 1/60 # degrees
dlon = 1/60 # degrees


# =============================================================================
# PROMPT
# =============================================================================
# PROMPT FOR GEOGRAPHIC REGION
ok = False
while not(ok):
    print("Enter limits (degrees): north south west east")
    print("latitudes (+90, -90)")
    print("longitudes (-180, +180), unless crossing prime meridian")
    dnorth, dsouth, dwest, deast = map( float, input().split(' ') )
    
    if ( (dnorth>90) or (dsouth<-90) or 
        (dwest>=deast) or (dsouth>=dnorth) or 
        (deast<-180) ):
        print("\t INVALID SELECTION")
    else:
        ok = True
print("\n")

# PROMPT FOR RESOLUTION
dlat_out = dlat
dlon_out = dlon

# PROMPT FOR OUTPUT FORMAT
#ok = False
#while not(ok):
#    print ("Output mode? (grid = 1, lat/lon/geoid_height = 2)")
#    iout = float(input())
#    
#    if ( (iout==1) or (iout==2) ):
#        ok = True
#    else: 
#        print("\t INVALID, enter 1 or 2")
#print("\n")
iout = 2; print("I'm only going for lat/lon/G_Height\n")


# =============================================================================
# OPEN INPUT AND OUTPUT FILES
# =============================================================================
fnul1  = f"{path_in}/{n_in}"
file_1 = FortranFile(fnul1, 'r')

fnul10 = f"{path_out}/{nout}"
file_10 = open(fnul10, "w+") # output file

print(f"Input Sequential Binary Data File : \n\t{n_in}\n")
print(f"Output Extracted Ascii Data File  : \n\t{nout}\n")
if (iout==1):
    print("ONE OUTPUT RECORD PER EXTRACTED PARALLEL BAND")
elif (iout==2):
    print("ONE OUTPUT RECORD PER EXTRACTED GRIDPOINT")
else:
    print("ERROR: IOUT MUST BE SET TO 1 OR 2")
print(line_5000)
print(line_5000)

     
# =============================================================================
# INITIALISE DIMENSIONS OF EXTRACTED GRID
# =============================================================================
north_i = round((90-dnorth) / dlat)
south_i = round((90-dsouth) / dlat) # included
west_j  = round((180+dwest) / dlon)
east_j  = round((180+deast) / dlon) # included

irow = south_i - north_i + 1
jcol = east_j - west_j + 1

if (west_j < 0):     west_j += ncols
if (west_j > ncols): west_j -= ncols

if (east_j < 0):     east_j += ncols
if (east_j > ncols): east_j -= ncols

north_m = 90 - north_i*dlat
south_m = 90 - south_i*dlat
west_m  = west_j*dlon - 180
east_m  = east_j*dlon - 180

print("Geometry of extracted grid:")
print(f" Latitude of northern-most points = {north_m} (Degrees)\n",
      f"Latitude of southern-most points = {south_m} (Degrees)\n",
      f"Longitude of western-most points = {west_m} (Degrees)\n",
      f"Longitude of eastern-most points = {east_m} (Degrees)\n\•n",
      f" Latitude spacing = {dlat_out*60} (Minutes)\n",
      f"Longitude spacing = {dlon_out*60} (Minutes)\n\n",
      f"{irow} Rows x {jcol} Columns of values output.")


# =============================================================================
# DIMENTIONS OK ?
# =============================================================================
if ((dnorth>90) or (dnorth<-90) or
   (dsouth>90) or (dsouth<-90) or
   (dsouth>dnorth) or (dwest>deast) or
   (deast<-180) or (jcol>ncols) ):
    print(line_5000)
    print("ERROR:CHECK INPUT DIMENSIONS")
    print(line_5000)


# =============================================================================
# EXTRACT GRID ONE ROW AT A TIME
# =============================================================================
print(line_5000)
print("Extracting values:\n")
        
flat = np.arange(north_m, south_m-dlat/2, -dlat_out)
flon = np.arange(west_m, east_m+dlon/2, dlon_out)   


# LOOP-READ THROUGH UNINTERESTING ROWS DATA 
for i in range(0, north_i):
    file_1.read_record(dtype=np.float32)


# LOOP-EXTRACT INTERESTING ROW DATA
for ii in range (0, irow):  
    grid = file_1.read_record(dtype=np.float32)
     
    temp = np.zeros(jcol)
    for jj in range (0, jcol):
        fetch_j = west_j +jj
        if (fetch_j>=ncols): fetch_j -= ncols # if you want fetch_j=21600
        print(f"fetch_i={north_i+ii}\t fetch_j={fetch_j}")
        temp[jj] = grid[fetch_j]
    
    if (iout==2):
        for j in range (0, jcol):
            file_10.write(f"{flat[ii]:.6f}\t{flon[j]:.6f}\t{temp[j]:.6f}\n")

stats()

file_1.close()
file_10.close()


# =============================================================================
# READ OUTPUT AND PRINT SAMPLE: TOP LEFT 10 PT X 10 PT SQUARE
# =============================================================================
if (iout==1):
    print(line_5000)
    print("Sample Output: Top Left 10pt x 10pt Square")
#    file_10 = open(fnul10, "r")    
    print("i didn't code this part")
    print(line_5000)

print(line_5000)
print("\t Normal Termination")
print(line_5000)


