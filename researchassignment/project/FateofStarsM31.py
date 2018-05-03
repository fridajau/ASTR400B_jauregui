import numpy as np
import astropy.units as u
import matplotlib
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from RingofStarsM31 import *

#--Compute the orbit of M31 at each pericenter and determine what percentage of stars will get
#transfered or pass through the milkyway <-- will need seperate data sets to tell the difference b/t particles from
#                                            the galaxies, maybe visual aid
#                                            if time allows M33 as well

#--Use lab7 to help rearrange the merger to be face-on and plot the density profile at snap shots
#--Create a histogram to see the radial distribution with repect to M31's galatic center of solar particles

def M31SunsOrbit(galaxy, Snap, end, n):
    #INPUT:  snap shots, integrator
    #OUTPUT: return the orbit of the M31 & MW system at deseried snap shots
    
    #--Compute the orbit of my ring particles with MW
    #--I have M31_Orbit.txt which contains the orbit of disk particles at all snap shots
    #--Use Orbits.py as a template
    #--Snaps shots to obtain (0.0, 3.87, 5.87, 6.2 & 10.0)
    #--file of the orbit of MW and the reduced M31 particles
    #--take steps from Orbits.py and hw5

    #--define an intinal array that will store a radii of particles
    Orbit = np.zeros((int(end/n)+1, 7))
    #--define the filename for the galaxy using hw5 as example
    for i in range(Snap,end+n,n):
        ilbl = '000' + str(i)
        ilbl = ilbl[-3:]
        filename = "%s_"%(galaxy) + ilbl + '.txt'
        #--COM M31 Sun-like stars from RingofStarsM31 at snapshot 0
        Disk_M31   = SolarParticles("M31",4, 2)
        radial_M31 = Disk_M31.RadialPos(-377,608,-284)
        COM   = CenterOfMass('VLowRes/'+filename, 2)
        #--store the COM pos and vel in the Orbit array, divide out the units
        Orbit[int(i/n),0] = float(COM.time/u.Myr)/1000
        Orbit[int(i/n),1] = radial_M31[0][:1]
        Orbit[int(i/n),2] = radial_M31[1][:1]
        Orbit[int(i/n),3] = radial_M31[2][:1]
        Orbit[int(i/n),4] = radial_M31[3][:1]
        Orbit[int(i/n),5] = radial_M31[4][:1]
        Orbit[int(i/n),6] = radial_M31[5][:1]
        #--check
        print(i)
            
    fileout = 'M31_Ring_0.txt'
    
    #--save the array Orbit to a file
    np.savetxt(fileout, Orbit, header='t, x, y, z, vx, vy, vz', comments='#',
                   fmt=['%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])
    
    return Orbit

M31_ring = M31SunsOrbit("M31", 0, 1, 1)


    #def CricVelM31():
        

