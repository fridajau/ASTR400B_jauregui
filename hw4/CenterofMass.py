"""
Frida Jauregui
Homework 4 started Feb 1, 2018
Center of Mass Position and Velocity
"""

import numpy as np
import numpy.ma as ma
import astropy.units as u
from ReadFile import Read

"""
Compute the center of mass posistion and velocity vectors of each galaxy at 
any given point.

Call the ReadFile function and the 3 text files.
"""
#---Class
class CenterOfMass:
   
    def __init__(self, filename, ptype):
        #---read in the file and particle type
        self.time, self.total, self.data = Read(filename)
            
        #---create an array to store indexes of particles of desired Ptype
        self.index = np.where(self.data['type'] == ptype)
    
        #---store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m'][self.index]
        
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]

        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]

#---Functions

#---Input:  x,y,z coords of the position OR velocity and mass
#---Output: returns the 3D coords of the COM(position OR velocity)

    def COMdefine(self, x,y,z, m):
        #---position
        pos = (np.sqrt((self.x)**2 + (self.y)**2 + (self.z)**2))*u.kpc
       
        #---velocity
        vel = (np.sqrt((self.vx)**2 + (self.vy)**2 + (self.vz)**2))*u.km/u.s

        #---total mass
        totmass = np.sum(self.m)*u.Msun*1e10
        
        #---while loop to return 3D position and total mass 
        while(self.total < 1):
            X = (pos*totmass)/totmass
            print(X)

        return pos, totmass

#---Input:  tolerance
#---Output: returns an array with the x,y,z coords of the COM position

    def COM_P(self, x,y,z, m):
        totmass = np.sum(self.m)*u.Msun*1e10
        
        #---COM position vector 
        XCOM = (self.x*totmass)/totmass
        YCOM = (self.y*totmass)/totmass
        ZCOM = (self.z*totmass)/totmass
        #---magnitude of the new position vector
        RCOM = np.linalg.norm(self.COMdefine(self, XCOM,YCOM,ZCOM, m))
        
        #---changing particle position vector 
        XNEW = (self.x - XCOM)
        YNEW = (self.y - YCOM)
        ZNEW = (self.z - ZCOM)

        RNEW = np.linalg.norm(self.COMdefine(self, XNEW,YNEW,ZNEW, m))
#######################################
        #---maximum 3D separation of the particles 
        RMAX = 0.5*(np.linalg.norm(self.x, self.y, self.z, m))

        #---while loop for RCOM - new COM > delta
        delta  = 10
        newCOM = 1000*u.kpc
        while(RCOM - newCOM > delta):
            RCOM2 = RCOM - newCOM
            print(RCOM2)/2

        return RCOM2

#---Input:  tolerance
#---Output: returns an array with the x,y,z coords of the COM position

    def COM_V(self, vx,vy,vz, m):
        totmass = np.sum(self.m)*u.Msun*1e10

        #---store velocities within 5kpc from the COM position
        vel      = (np.sqrt((self.vx)**2 + (self.vy)**2 + (self.vz)**2))*u.km/u.s
        newindex = ma.masked_where(vel > 15, vel)
        newvel   = vel[newindex]

        VCOM = self.COMdefine(self, vel, m)

        return VCOM
    
#---prints COM position and velocity of MW disk properties
MWCOM = CenterOfMass("MW_000.txt", 2)
#MW_P  = MWCOM.COM_P(0,0,0,0)
#print("MW COM Disk:", MWCOM)

#---prints vaule of disk mass of galaxy and the 3d position vector              
MW_mass = MWCOM.COMdefine(0,0,0,0)
print("MW Disk Mass:", MW_mass)

M31COM = CenterOfMass("M31_000.txt", 2)
M31_mass = M31COM.COMdefine(0,0,0,0)
print("M31 Disk Mass:", M31_mass)

M33COM = CenterOfMass("M33_000.txt", 2)
M33_mass = M33COM.COMdefine(0,0,0,0)
print("M33 Disk Mass:", M33_mass)

