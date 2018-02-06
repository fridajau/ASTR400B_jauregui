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

    def total_mass(self):
        return np.sum(self.m)*u.Msun*1e10

#---Input:  x,y,z coords of the position OR velocity and mass
#---Output: returns the 3D coords of the COM(position OR velocity)

    def COMdefine(self, x,y,z, m):
        #---position and total mass
        #pos3d = (np.sqrt(self.x**2 + self.y**2 + self.z**2))*u.kpc
        totmass = np.sum(m)
        #---define COM of x,y,z
        XCOM = np.sum(x*m)/(totmass)
        YCOM = np.sum(y*m)/(totmass)
        ZCOM = np.sum(z*m)/(totmass)
    
        return XCOM,YCOM,ZCOM

#---Input:  tolerance (small like 1)
#---Output: returns an array with the x,y,z coords of the COM position

    def COM_P(self, delta):
        #---COM position vector estimate
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        
        #---magnitude of the desired vector
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)
        
        #---changing particle position vector to the COM frame 
        XNEW = (self.x - XCOM)
        YNEW = (self.y - YCOM)
        ZNEW = (self.z - ZCOM)
        #---mangitude of the new position vectors
        RNEW = np.sqrt(XNEW**2 + YNEW**2 + ZNEW**2)
        
        #---maximum 3D separation of the particles from the COM position
        RMAX = np.max(RNEW)/2

        #---while loop for RCOM - new COM > delta
        diff = 1000
        while(diff) > delta:
            #---index particles 
            newindex = np.where(RMAX > RNEW)
            x2 = self.x[newindex]
            y2 = self.y[newindex]
            z2 = self.z[newindex]
            m2 = self.m[newindex]
            
            #---redefine for 1st guess
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2,y2,z2, m2)
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)
            
            #---redefine for new change of position vector
            ncx = (x2 - XCOM2)
            ncy = (y2 - YCOM2)
            ncz = (z2 - ZCOM2)
            RNEW2 = np.sqrt(ncx**2 + ncy**2 + ncz**2)
            diff = np.abs(RCOM2 - RCOM)
            #---get new max seperation
            RMAX2 = np.max(RNEW2)/2

        return XCOM2,YCOM2,ZCOM2

#---Input:  vx,vy,vz and mass
#---Output: returns an array with the x,y,z coords of the COM position

    def COM_V(self, vx,vy,vz):
        #---store velocities
        VXCOM, VYCOM, VZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        VRCOM = np.sqrt(VXCOM**2 + VYCOM**2 + VZCOM**2)

        #---store velocities within 15kpc from the COM position
        vindex = np.where(VRCOM > 15)
        vxnew  = vx[vindex]
        vynew  = vy[vindex]
        vznew  = vz[vindex]
        VXCOM2, VYCOM2, VZCOM2 = (vxnew,vynew,vznew)
        
        return VXCOM2,VYCOM2,VZCOM2 
    

#---prints vaule for x,y,z
MWCOM = CenterOfMass("MW_000.txt", 2)
MW_mass = MWCOM.COMdefine(MWCOM.x, MWCOM.y, MWCOM.z, MWCOM.m)
print("MW Disk Mass:", MW_mass)

M31COM = CenterOfMass("M31_000.txt", 2)
M31_mass = M31COM.COMdefine(M31COM.x, M31COM.y, M31COM.z ,M31COM.m)
print("M31 Disk Mass:", M31_mass)

M33COM = CenterOfMass("M33_000.txt", 2)
M33_mass = M33COM.COMdefine(M33COM.x, M33COM.y, M33COM.z, M33COM.m)
print("M33 Disk Mass:", M33_mass)

#---Testing your code
MW_pos = MWCOM.COM_P(1)
print("COM Disk position for MW:", MW_pos)

M31_pos = M31COM.COM_P(1)
print("COM Disk position for M31:", M31_pos)

M33_pos = M33COM.COM_P(1)
print("COM Disk position for M33:", M33_pos)


