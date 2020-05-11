# -*- coding: utf-8 -*-
"""
Created on Sun May 10 01:33:46 2020

@author: Xavier

Fortran code for gridget xmin

c-----------------------------------------------------------------------
c     ORIGINAL PROGRAM:                           SIMON HOLMES, JUL 2007
c     MODIFIED FOR CORNER-CELL REGISTRATION       SIMON HOLMES, MAY 2008
c     MUTLI-OUTPUT OPTION                         SIMON HOLMES, MAY 2008
#     TRANSLATION INTO PYTHON              XAVIER DE LABRIOLLE, MAY 2020
c-----------------------------------------------------------------------

"""
import numpy as np
from ast import literal_eval
import struct
from scipy.io import FortranFile

term_width = 50
line_5000 = "-" * term_width

def nint(x): # OR THE round(à) function
#    i = int(x)
#    if (abs(x-i)>=0.5): i +=1
#    return i
    return round(x)

def stats():
    pass

# =============================================================================
# =============================================================================
# # MAIN
# =============================================================================
# =============================================================================


# ---- line 78
# =============================================================================
# INPUT PARAMETERS
# =============================================================================
path_in  = './Fortran grid'
#n_in     = 'Und_min1x1_egm2008_isw=82_WGS84_TideFree_SE'
#n_in     = 'RAW_binary_test'
n_in     = 'RAW_binary_test2.txt'
path_out = '.'
nout     = 'OUTPUT.DAT'

dlat = 1/60 # degrees
dlon = 1/60 # degrees

# =============================================================================
# NON-INPUT PARAMETERS
# =============================================================================
nrows = 10801 # EGM2008 grid size
ncols = 21600

scrd  = np.zeros((ncols,2))
statd = np.zeros(22)
flat  = np.zeros((nrows))
flon  = np.zeros(ncols)
grid  = np.zeros(ncols)
temp  = np.zeros(ncols)

print(line_5000)
print("\t Execution")
print(line_5000)

'''
# ---- 109
# =============================================================================
# PROMPT OR GEOGRAPHIC REGION
# =============================================================================
ok = False
while not(ok):
    print("Enter limits(degrees): north south west east")
    print("longitudes 0-360, unless crossing prime meridian")
    dnorth, dsouth, dwest, deast = map( float, input().split(' ') )
    
    if ( (dnorth>90) or (dsouth<-90) or 
        (dwest>=deast) or (dsouth>=dnorth) or 
        (deast<0) ):
        print("\t INVALID SELECTION")
    else:
        ok = True
    print("\n")

# =============================================================================
# PROMPT FOR RESOLUTION
# =============================================================================
print("Not yet asking for resolution")

# =============================================================================
# PROMPT FOR OUTPUT FORMAT
# =============================================================================

ok = False
while not(ok):
    print ("Output mode? (grid = 1, lat/lon/geoid_height = 2)")
    iout = float(input())
    
    if ( (iout==1) or (iout==2) ):
        ok = True
    else: 
        print("\t INVALID, enter 1 or 2")
    print("\n")

iout = 2; print("I'm only going for lat/lon/G_Height anyway bitches\n")

'''
# ---- line 152
# =============================================================================
# OPEN INPUT AND OUTPUT FILES
# =============================================================================
fnul1  = f"{path_in}/{n_in}"
file_1  = open(fnul1, "rb") # text file with binary coordinates
#field = np.fromfile(file_1, dtype='float32')
#    file_1.close(); file_1  = open(fnul1, "rb")

file_1 = FortranFile(fnul1, 'r')








'''
header = file_1.read(6)
arr = []
for i in range (0, 21600):
    arr.append( struct.unpack('f', file_1.read(4))[0] )


#    print(arr[-1])

txt = file_1.read(50)
print(txt)


MAX = len(arr)

ii = 0
#while (ii<MAX):
aa = arr[ii]
while (arr[ii]==aa): ii += 1
print(f"all equal to {aa} until ii={ii}")


iii=21600

points = 25
for jj in range (iii, iii+points):
    print(f"point {jj} is {arr[jj]}")




for y in range (0, 10):
    ct = -1
    file_1.close()
    file_1  = open(fnul1, "rb")
    
    header = file_1.read(6)
    ct +=6
    
    arr = []
    for i in range (0, ncols):
        arr.append( struct.unpack('f', file_1.read(4))[0] );
        ct+=4
    
    file_1.read(y)
    ct += y
    arr.append( struct.unpack('f', file_1.read(4))[0] ); ct +=4
    
    print(f"ct={ct}\ty={y}\tdata={arr[-1]}")
    input()


for y in range(14, 15):
    ct = -1
    file_1.close()
    file_1  = open(fnul1, "rb")
    
    header = file_1.read(6)
    ct +=6
    
    file_1.read(ncols*4)
    ct += ncols*4
    
    file_1.read(y)
    ct += y
    
    arr = []
    for i in range (0, ncols):
        arr.append( struct.unpack('f', file_1.read(4))[0] )
        file_1.read(2)
        ct+=6
        

#    arr.append( struct.unpack('f', file_1.read(4))[0] ); ct +=4
    
    print(f"ct={ct}\ty={y}\tdata={arr[-1]}")
    input()
'''










'''
fnul10 = f"{path_out}/{nout}"
file_10 = open(fnul10, "w+") # output file

print(f"Input Sequential Binary Data File : \n\t{n_in}")
print(f"Output Extracted Ascii Data File  : \n\t{nout}")
print("\n")
if (iout==1):
    print("ONE OUTPUT RECORD PER EXTRACTED PARALLEL BAND")
elif (iout==2):
    print("ONE OUTPUT RECORD PER EXTRACTED GRIDPOINT")
else:
    print("ERROR: IOUT MUST BE SET TO 1 OR 2")
print(line_5000)
print(line_5000)


# ---- line 189     
# =============================================================================
# INITIALISE DIMENSIONS OF EXTRACTED GRID
# =============================================================================
d_8 = 1e-8
iglob = 1 - (int(((deast-dwest)/360) + d_8))

in_ = round( ((90-dnorth)/dlat) - d_8) + 1
is_ = round( ((90-dsouth)/dlat) - d_8) + 1
jw  = round(dwest/dlon - d_8) + 1
je  = round(deast/dlon + d_8) + iglob
irow = is_ - in_ + 1
jcol = je - jw + 1

jwl = jw
if (jwl < 1):     jwl += ncols
if (jwl > ncols): jwl -= ncols

jel = je
if (jel < 1):     jel += ncols
if (jel > ncols): jel -= ncols

tlat = 90 - (in_ - 1)*dlat
wlon =      (jwl - 1)*dlon
slat = 90 - (is_ - 1)*dlat
elon =      (jel - 1)*dlon

print("Geopetry of extracted grid:")

print(f" Latitude of northern-most points = {tlat} (Degrees)\n",
      f"Latitude of southern-most points = {slat} (Degrees)\n",
      f"Longitude of western-most points = {wlon} (Degrees)\n",
      f"Longitude of eastern-most points = {elon} (Degrees)\n",
      f"                Latitude spacing = {dlat*60} (Minutes)\n",
      f"               Longitude spacing = {dlon*60} (Minutes)\n\n",
      f"{irow} Rows x {jcol} Columns of values output.")

dn = 90 - (in_ - 1)*dlat # for stats sub only
dw =      (jw  - 1)*dlon

# =============================================================================
# DIMENTIONS OK ?
# =============================================================================
if ((dnorth>90) or (dnorth<-90) or
   (dsouth>90) or (dsouth<-90) or
   (dsouth>dnorth) or (dwest>deast) or
   (deast<0) or (jcol>ncols) ):
    print(line_5000)
    print("ERROR:CHECK INPUT DIMENSIONS")
    print(line_5000)


# ---- line 248
# =============================================================================
# EXTRACT GRID ONE ROW AT A TIME; COMPUTE STATS FOR EXTRACTED GRID
# =============================================================================
if (iout==1):
    print("please try again with iout = 2 thank you")
# ----
if (iout==2):
    for i in range (0, nrows):
        flat[i] = 90 - (i)*dlat
# ----
for i in range(0, in_ - 1): # reaching the point we are interested in
    file_1.read(1)

print(line_5000)
print("Statistics of Extracted values:\n")

ii = 0
for i in range (in_, is_+1):
    file_1.read(1)
    ii += 1
    jj = 0
    for j in range (jw, je+1):
        jl = j
        if (jl<1):     jl += ncols
        if (jl>ncols): jl -= ncols
        jj += 1
        temp[jj] = grid[jl]
        scrd[jj, 0] = grid[jl]
        scrd[jj, 1] = 1
        flon[jj] = (jl-1)*dlon
    
    if (iout==1): 
        file_10.write(str(temp[j]))
    
    if (iout==2):
        for j in range (0, jj-1):
            file_10.write(f"{flat[i]}\t{flon[j]}\t{temp[j]}\n")
    
    stats()

file_10.close()


# ---- line 301
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

'''
