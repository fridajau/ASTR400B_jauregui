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
        pos3d = (np.sqrt(self.x**2 + self.y**2 + self.z**2))*u.kpc
       
        #---total mass
        totmass = np.sum(self.m)*u.Msun*1e10
        
        #---define COM of x,y,z
        XCOM = np.sum(self.x*m)/(totmass)
        YCOM = np.sum(self.y*m)/(totmass)
        ZCOM = np.sum(self.z*m)/(totmass)
        pos = (XCOM,YCOM,ZCOM)

        return pos, totmass

#---Input:  tolerance (small like 1)
#---Output: returns an array with the x,y,z coords of the COM position

    def COM_P(self, delta):
        #---COM position vector estimate
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        
        #---magnitude of the desired vector
        RCOM = np.linalg.norm(self.COMdefine(XCOM,YCOM,ZCOM, m))
        
        #---changing particle position vector to the COM frame 
        XNEW = (self.x - XCOM)
        YNEW = (self.y - YCOM)
        ZNEW = (self.z - ZCOM)
        #---mangitude of the new position vectors
        RNEW = np.linalg.norm(self.COMdefine(XNEW,YNEW,ZNEW, m))
        
        #---maximum 3D separation of the particles from the COM position
        RMAX = np.max(RNEW)/2

        #---while loop for RCOM - new COM > delta
        RCOM2 = 1000*u.kpc
        while(RCOM2) > delta:
            #---index particles 
            newindex = np.where(RMAX > RNEW)
            x2 = self.x[newindex]
            y2 = self.y[newindex]
            z2 = self.z[newindex]
            m2 = self.m[newindex]
            #---redefine for 1st guess
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(self.x, self.y, self.z, self.m)
            RCOM2 = np.linalg.norm(self.COMdefine(XCOM2,YCOM2,ZCOM2, m))
            
            #---new change of position vector
            ncx = (x2 - XCOM2)
            ncy = (y2 - YCOM2)
            ncz = (z2 - ZCOM2)

            
            #RCOM2 = RCOM - COM
            #print(RCOM2)/2

        return RCOM2

#---Input:  vx,vy,vz and mass
#---Output: returns an array with the x,y,z coords of the COM position

    def COM_V(self, vx,vy,vz, m):
        totmass = np.sum(self.m)*u.Msun*1e10

        #---store velocities within 5kpc from the COM position
        vel      = (np.sqrt((self.vx)**2 + (self.vy)**2 + (self.vz)**2))*u.km/u.s
        newindex = ma.masked_where(vel > 15, vel)
        newvel   = vel[newindex]

        VCOM = self.COMdefine(newvel, m)

        return VCOM





    

#---prints vaule of disk mass of galaxy
MWCOM = CenterOfMass("MW_000.txt", 2)
MW_mass = MWCOM.COMdefine(0,0,0,0)
print("MW Disk Mass:", MW_mass)

M31COM = CenterOfMass("M31_000.txt", 2)
M31_mass = M31COM.COMdefine(0,0,0,0)
print("M31 Disk Mass:", M31_mass)

M33COM = CenterOfMass("M33_000.txt", 2)
M33_mass = M33COM.COMdefine(0,0,0,0)
print("M33 Disk Mass:", M33_mass)

