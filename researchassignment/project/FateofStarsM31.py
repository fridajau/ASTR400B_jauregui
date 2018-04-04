import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterofMass import CenterOfMass
from Orbits import OrbitCOM



#---class
#---index disk particles from data to work only with disk particles
class SolarParticles:
    def __init__(self, filename, ptype):
        self.time, self.total, self.data = Read(filename)
        #--storing data for disk particles
        self.index = np.where(self.data['type'] == ptype)
        
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]

#--Index data to return vaules of M31 within a radial distance of 7-9 kpc.
#--Contain  particles that are within half the scale hieght of M31
#--Get COM of the orbits of the galaxies and change the particle pos/vel vector to the COM frame for
#each snap shot <--obtain files from this
#--Compute the circular velocity of M31 at each pericenter and determine what percentage of stars will get
#transfered or pass through the milkyway <-- will need seperate data sets to tell the difference b/t particles from
#                                            the galaxies, maybe visual aid
#                                            if time allows M33 as well

#--Use lab7 to help rearrange the merger to be face-on and plot the density profile at snap shots
#--Create a histogram to see the radial distribution with repect to M31's galatic center of solar particles

    def RadialPos(self, x,y,z, vx,vy,vz):
        #INPUT:   x,y,z pos vx,vy,vz vel
        #RETURNS: 3d coord of pcom and vcom of particals within 7-9 kpc

        #--Use Orbits file to obtain x,y,z vx,xy,vz COM (0.0, 3.87, 5.87, 6.2, 10.0)
        #--genfrom txt the Orbits of MW and M31
        #--index M31 and return an array with the desired particles
        MW  = np.genfromtxt('MW_Orbit.txt', dtype = float, names=True)
        M31 = np.genfromtxt('M31_Orbit.txt', dtype = float, names=True)

        #--select particles from m31
        tm31 = M31['t']
        xm31 = M31['x']
        ym31 = M31['y']
        zm31 = M31['z']

        vxm31 = M31['vx']
        vym31 = M31['vy']
        vzm31 = M31['vz']

        #--radial position from galatic center
        R = np.sqrt(xm31**2 + ym31**2)
        
        Rindex = np.where((R > 7) & (R < 9))
        x2m31 = xm31[Rindex]
        y2m31 = ym31[Rindex]
        z2m31 = zm31[Rindex]

        vx2m31 = vxm31[Rindex]
        vy2m31 = vym31[Rindex]
        vz2m31 = vzm31[Rindex]

        
        
