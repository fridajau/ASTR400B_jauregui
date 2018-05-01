import numpy as np
import astropy.units as u
import matplotlib
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterofMass import CenterOfMass

#---class
#---index disk particles from data to work only with disk particles
class SolarParticles:
    def __init__(self, galaxy, snap, ptype):
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        # create filenames
        self.filename='%s_'%(galaxy) + ilbl + '.txt'
        self.time, self.total, self.data = Read('VLowRes/'+self.filename)
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

    def RadialIndex(self, x,y,z):
        #INPUT:  x,y,z postions
        #OUTPUT: the ring of particles

        #--use CenterOfMass to obtain the COM of x,y,z vx,xy,vz of M31
        #--snap shots to look out for (0.0, 3.87, 5.87, 6.2, 10.0)
        #--obtain the index ring of particles to use later on to define at wanted snap shots
        #--follow this function for radial positions

        #--creating a COM obeject for M31 using Disk Particles from CenterOfMass
        COM_M31 = CenterOfMass('VLowRes/'+self.filename,2)
        M31_pos = COM_M31.COM_P(1.0,4.0)
        M31_vel = COM_M31.COM_V(M31_pos[0],M31_pos[1],M31_pos[2])
        #--x,y,z pos & vel of m31 rel to COM
        COMX = COM_M31.x - float(M31_pos[0]/u.kpc)
        COMY = COM_M31.y - float(M31_pos[1]/u.kpc)
        COMZ = COM_M31.z - float(M31_pos[2]/u.kpc)
        
        COMVX = COM_M31.vx - float(M31_vel[0]/(u.km/u.s))
        COMVY = COM_M31.vy - float(M31_vel[1]/(u.km/u.s))
        COMVZ = COM_M31.vz - float(M31_vel[2]/(u.km/u.s))
        #--3d radial position from galatic center
        RadPos = np.sqrt(COMX**2 + COMY**2 + COMZ**2)
        #--index the vaules within 7-9 kpc 
        Rindex = np.where((RadPos > 7) & (RadPos < 9) & (COMZ > -18.5) & (COMZ < 18.5))
    
        return Rindex
    
    def RadialPos(self, x,y,z):
        #INPUT: x,y,z positions, center of galaxy=0  
        #RETURNS: 3d coord of pcom or vcom of particals within 7-9 kpc a snap shot 0

        #--index M31 vaules and return an array with the desired particles

        COM_M31 = CenterOfMass('VLowRes/'+self.filename,2)
        M31_pos = COM_M31.COM_P(1.0,4.0)
        M31_vel = COM_M31.COM_V(M31_pos[0],M31_pos[1],M31_pos[2])
       
        COMX = COM_M31.x - float(M31_pos[0]/u.kpc)
        COMY = COM_M31.y - float(M31_pos[1]/u.kpc)
        COMZ = COM_M31.z - float(M31_pos[2]/u.kpc)
        
        COMVX = COM_M31.vx - float(M31_vel[0]/(u.km/u.s))
        COMVY = COM_M31.vy - float(M31_vel[1]/(u.km/u.s))
        COMVZ = COM_M31.vz - float(M31_vel[2]/(u.km/u.s))

        #--3d radial position from galatic center
        RadPos = np.sqrt(COMX**2 + COMY**2 + COMZ**2)
        #--use Rindex function to always return the ring of particles
        Rindex = self.RadialIndex(COMX,COMY,COMZ)
        nM31x = COMX[Rindex]
        nM31y = COMY[Rindex]
        nM31z = COMZ[Rindex]

        nM31vx = COMVX[Rindex]
        nM31vy = COMVY[Rindex]
        nM31vz = COMVZ[Rindex]

        return nM31x, nM31y, nM31z, nM31vx, nM31vy, nM31vz

#--testing code
Disk_M31 = SolarParticles("M31", 3, 2)
radial_ring = Disk_M31.RadialPos(-377,608,-284)
#radial_index = Disk_M31.RadialIndex(-377,608,-284)

testx  = np.around(radial_ring[0],3)
testy  = np.around(radial_ring[1],3)
testz  = np.around(radial_ring[2],3)
testvx = np.around(radial_ring[3],3)
testvy = np.around(radial_ring[4],3)
testvz = np.around(radial_ring[5],3)

fewpart = testx[:3]
print(fewpart)

rad = np.average(np.sqrt(testx**2 + testy**2 + testz**2))
vel = np.average(np.sqrt(testvx**2 + testvy**2 + testvz**2))

print("xpos at snap 003", testx)
#print("ave pos", np.around(rad,3))
#print("ave vel", np.around(vel,3))
