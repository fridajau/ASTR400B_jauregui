import numpy as np
import astropy.units as u
import matplotlib
from matplotlib.colors import LogNorm
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

    def RadialPos(self, x,y,z):
        #INPUT: x,y,z positions, center of galaxy=0  
        #RETURNS: 3d coord of pcom or vcom of particals within 7-9 kpc

        #--Use Orbits txt to obtain the COM of x,y,z vx,xy,vz of M31
        #--snap shots to look out for (0.0, 3.87, 5.87, 6.2, 10.0)
        #--index M31 vaules and return an array with the desired particles

        #--create a COM obeject for M31 using Disk Particles from CenterOfMass
        COM_M31 = CenterOfMass("M31_000.txt",2)
        M31_pos = COM_M31.COM_P(1.0,4.0)
        M31_vel = COM_M31.COM_V(M31_pos[0],M31_pos[1],M31_pos[2])

        #--x,y,z pos & vel of m31 rel to COM
        COMX = COM_M31.x - float(M31_pos[0]/u.kpc)
        COMY = COM_M31.y - float(M31_pos[1]/u.kpc)
        COMZ = COM_M31.z - float(M31_pos[2]/u.kpc)

        #--3d radial position from galatic center
        RadPos = np.sqrt(COMX**2 + COMY**2)

        #--index the vaules within 7-9 kpc
        Rindex = np.where((RadPos > 7) & (RadPos < 9))
        ###Zindex = np.where((zm31 > -1.5) & (zm31 < 1.5))
        nM31x = COMX[Rindex]
        nM31y = COMY[Rindex]
        nM31z = COMZ[Rindex]
        ###z3m31 = z2m31[Zindex]
                
        #--particals within the new radius    
        RadPos2 = np.sqrt(nM31x**2 + nM31y**2)
        return RadPos2
    
#--work with disk particles of M31
Disk_M31   = SolarParticles("M31_000.txt", 2)
#--take x,y,z to be the relative position of M31 to MW from CenterOfMass
SunCandM31 = self.RadialPos(-377, 608, -284)
    
    #--Compute the orbit of my ring particles with MW
    #--I have MW_Orbit.txt which contains the orbit of disk particles of all snap shots
    #--Use Orbits.py as a template
    #--Snaps shots to obtain (0.0, 3.87, 5.87, 6.2 & 10.0)
    
    def M31SunsOrbit(self, DiskPart, Snap):
        #INPUT:  the particles I'm going to orbit, snap shots
        #OUTPUT: return the orbit of the M31 & MW system at deseried snap shots
        
        #--file of the orbit of MW and the reduced M31 particles
        #--take steps from Orbits.py, if/when statement so not rerun code over and over 
        fileout = "Orbit_{}.txt".format(galaxy)

    



    

        
