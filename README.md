# tolosat_grav_payload
This is the github repository for the software tools developped by the gravimetry payload of the TOLOSAT CubeSat. 
For more information on TOLOSAT, follow us on twitter at: https://twitter.com/tolosat1

# ===============================
# ======== Grav Harm 3.0 ========
# ===============================


# ==== About this tool ====
This tool has been developped to strengthen our claims, here at the TOLOSAT Gravimetry Payload, for phases A and B of this CubeSat mission. I have late-night-baptized this tool:
	"Grav Harm 3.0" 
for "Gravitational Harmonics 3.0". 
It is the third attemp at writing this tool from the ground up. 
Any other name ideas are gladly welcome.

This tool was written in python 3.6, since 2019, by Xavier de Lariolle, [add your names if you've contributed, folks]. 


# ==== What does it do? ====
The purpose of this tool is to solve for the Stokes "spherical harmonic coefficients" to the geopotential or Earth. 
The calculation process is as follows: 
	1. Obtain the positions to the trajectory of a satellite
	2. Derivate the position to find the acceleration
	3. Solve for the coefficients, by various methods, of the geopotential
	4. Calculate the height of the Geoid undulation over the surface of the Earth
	5. Compare the calculated model with previous works in scientiic litterature
	6. Generate specific satellite ephemeris to answer the mission questions and confirm the mission requirements	


# ==== Before you start ====
In order to work properly, this tool requires an additional directory of files called "data", that can be downloaded at:
https://drive.google.com/drive/folders/19jqZdLf2ABAerFzBV8e-6N0IGVGUOotc
This directory contains coefficient files and simulated satellite trajectory files. It must be placed in the git repository, sams directory as "source", which contains all the tool's precious scripts. If you cannot access this Google Drive link, you aren't supposed to anyway (for now). If you really want access to it, send an email to any of the contributors.

The simulated satellite ehpemeris have been rendered using NASA's free open source GMAT software.
More on GMAT at: https://software.nasa.gov/software/GSC-17177-1

To download git for Windows, use this guide: https://www.computerhope.com/issues/ch001927.htm


# ==== Getting started ====
I suggest using anaconda to download libraries and spyder3 as an IDE.

The neded libraries are: numpy, scipy, matplotlib, cartopy, os, time, math, [add if you use more]. 
This tool contains many scripts, each regrouping the functions related to a specific task. The name of the task is reflected in the name of the script. 

It all happens in the script "GH_user". 
import the ephemeris and coefficient files you desire, 
set the maximum degree wanted for the calculations, 
set the plotting parameters, set the saving parameters.

The user must have some knowledge of python to write scripts, but many basic functions hae aleady been written and can be called from other scripts.


# ==== Just the Geoid please ====
IF ALL YOU WANT TO DO IS CALCULATE THE GEOID ELEVATION AT LAT/LONG COORDINATES, I'M NOT DONE CODING IT AND THE FUNCTIONS DO NOT MAKE SENSE YET. GIVE ME SOME TIME, THANKS. 


# ==== Sources and material ====
It is important to understand the concepts behing the geopotential, the geoid, the refference ellipsoid, spherical harmonics, Stokes coefficients. 

The main sources for the mathematics involved in this code are:

	"Definition of Functionals of the Geopotential and Their Calculation from Spherical Harmonic Models"
	by Franz Barthelmes
	for ICGEM
	Revision: jan. 2013
	Will be refered to as "GFZ".
	
	"How to Compute Geoid Undulations (Geoid Height Relative to a Given Reference Ellipsoid) from Spherical Harmonic Coefficients for Satellite Altimetry Applications"
	by Martin Losch and Verena Seufer
	dates back to 2003
	Heavily relies on the textbooks" Heiskanen & Moritz (1967)", and "Torge (1991)".
	Will be refered to as "The Geoid Cook Book"


# ==== Some information for coders and nerds ====  
To calculate geopotential values and to describe satellite trajectories, this python tool uses the spherical coordinates system :
 
	r = radius in km ; [0, inf()]
	theta = inclination from the z axis in radians; [0, pi]
	phi = rotation around the z axis in radians; [0, 2*pi]
 
They must be adapted to lat/lon coordinates by: 
	Lat = (pi/2 - theta) * 180/pi
	Long = (phi - pi) * 180/pi
 
Their equivalent in the cartesian coordinates is: 
	x = r * sin(theta) * cos(phi)
	y = r * sin(theta) * sin(phi)
	z = r * cos(theta)
	
	This is function: x,y,z = GH_convert.sph2cart(r,theta,phi)


# ==== Credits ====
Xavier C. de Labriolle
Antoine Bendimerad
Cedric Belmant

