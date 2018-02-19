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
        
        #--create a COM object using disk particles
        COM   = CenterOfMass(filename, 2)
        GCOMP = COM.COM_P(1.0, 2.0)
        GCOMV = COM.COM_V(GCOMP[0],GCOMP[1],GCOMP[2])


        #--first column in Orbit, store the time in Gyr(divide by 1000)
        #--row index of the Orbit array given as int(i)/n)
        #--quantities cant be stored in arrays divide out the units
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
        
    #--rerun code for MW,M31,M33    
    fileout = 'MW_lifetime.txt'
    
    #--save the array Orbit to a file
    np.savetxt(fileout, Orbit, header='t, x, y, z, vx, vy, vz', comments='#',
                   fmt=['%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])

    return Orbit

#--check code
MW = OrbitCOM('MW',0,800,5)

"""    
#---plots maybe in seperate python code
#--compute the time, comp and comv for each galaxy from snap 0 to 800
#--generate text files
MW = np.genfromtxt('MW_lifetime.txt', dtype = float, names=True)

###POS
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

plt.xlabel('Time')
plt.ylabel('X position')
plt.plot(MW['t'], MW['x'], label='MW', color='red')
plt.plot(M33['t'], M33['x'], label='M33', color='green')
plt.plot(M31['t'], M31['x'], label='M31', color='blue')

plt.show()
"""
