##################################
#MODIFIED ONLY FILEOUT
##################################
"""
Frida Jauregui
Homework 6 started Feb. 16, 2018

store the seperation and relative vel of MW, M31 & M33 over 
the entire simulation and plot the corresponding orbits

compute orbits up to snapshot 800(~12Gyr) of intervals of n=5

"""

import numpy as np
import astropy.units as u
from ReadFile import Read
from CenterofMass import CenterOfMass
import matplotlib.pyplot as plt
import os, glob 
#---create a symbolic link to the dir
#ln -s /home/astr400b/VLowRes /home/fjauregui


def OrbitCOM(galaxy, start, end, n):
    #Input: galaxy name
    #-------number of the 1st snapshot to be read in
    #-------the last snapshot to be read in
    #-------integer indicating the intervals to return COM
    #Output:file of time & COM pos and vel vectors of a given galaxy
    
    
    #--define the filename for the file that will store the orbit
    fileout = "Orbit_{}.txt".format(galaxy)

    #--define an array with 7 columns that will store
    #--x,y,z vx,vy,vz of the COM of the galaxy at each snapshot
    Orbit   = np.zeros((int(end/n)+1, 7))
    
    #--set delta & VolDec to be used CenterOfMass object for disk particles
    #--as in hw5 disk Particles afford the best centroiding
    #--M33 these vaules need to be small: delta=0.5 and VolDec=4
    delta  = 1
    Voldec = 2

    #--for loop from start to end+1 in intervals of n
    for i in np.arange(start, end+n, n):
        #--define the filename for the galaxy using hw5 as example
        ilbl = '000' + str(i)
        ilbl = ilbl[-3:]
        filename = "%s_"%(galaxy) + ilbl + '.txt'
        
        #--must add VLowRes txt files for 'VLowRes/' and filename
        #os.chdir("/home/fjauregui/VLowRes/")
        #filename2 = glob.glob("*.txt")
        
        #--create a COM object using disk particles
        COM   = CenterOfMass('VLowRes/'+filename, 2)
        GCOMP = COM.COM_P(1.0, 2.0)
        GCOMV = COM.COM_V(GCOMP[0],GCOMP[1],GCOMP[2])


        #--first column in Orbit, store the time in Gyr(divide by 1000)
        #--row index of the Orbit array given as int(i)/n)
        #--quantities cant be stored in arrays divide out the units
        #time       = float(float(COM.time/u.Myr)/1000)
        #store_time = np.insert(Orbit,1,time, axis=0)
        Orbit[int(i/n),0] = float(COM.time/u.Myr)/1000
        
        #--store the COM pos and vel in the Orbit array, divide out the units
        Orbit[int(i/n),1] = float(GCOMP[0]/u.kpc)
        Orbit[int(i/n),2] = float(GCOMP[1]/u.kpc)
        Orbit[int(i/n),3] = float(GCOMP[2]/u.kpc)


        Orbit[int(i/n),4]= float(GCOMV[0]/(u.km/u.s))
        Orbit[int(i/n),5]= float(GCOMV[1]/(u.km/u.s))
        Orbit[int(i/n),6]= float(GCOMV[2]/(u.km/u.s))

        #Orbit[int(i/n),0] = np.array([time, xcom,ycom,zcom, vxcom,vycom,vzcom])

        #--print the counter for the loop to the screen to kno where the code at
        print(i)
        
    #--rerun code for MW,M31,M33 with changed delta and VolDec vaules   
    fileout = 'M31_Orbit.txt'
    
    #--save the array Orbit to a file
    np.savetxt(fileout, Orbit, header='t, x, y, z, vx, vy, vz', comments='#',
                   fmt=['%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])

    return Orbit

#--check code
#MW  = OrbitCOM('MW',0,800,1)
#M31 = OrbitCOM('M31',0,800,1)





