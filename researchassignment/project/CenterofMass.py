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

    def COMdefine(self, a,b,c, m):
        #---define COM of x,y,z
        Acom = np.sum(a*m)/np.sum(m)
        Bcom = np.sum(b*m)/np.sum(m)
        Ccom = np.sum(c*m)/np.sum(m)
    
        return Acom,Bcom,Ccom

#---Input:  tolerance (small like 1)
#---Output: returns an array with the x,y,z coords of the COM position

    def COM_P(self, delta, VolDec):
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
        RMAX = max(RNEW)/VolDec

        #---while loop for RCOM - new COM > delta
        diff = 1000
        while(diff > delta):
            #---index particles 
            newindex = np.where(RNEW < RMAX)
            x2 = self.x[newindex]
            y2 = self.y[newindex]
            z2 = self.z[newindex]
            m2 = self.m[newindex]
            
            #---redefine for 1st guess
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2,y2,z2, m2)
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            #---determine the diff b/t the previous COM
            diff = np.abs(RCOM - RCOM2)
            #---reduce the volume by a factor of 2 again
            RMAX = RMAX/VolDec
            
            #---redefine for new change of position vector
            XNEW = (self.x - XCOM2)
            YNEW = (self.y - YCOM2)
            ZNEW = (self.z - ZCOM2)
            RNEW = np.sqrt(XNEW**2 + YNEW**2 + ZNEW**2)

            #---set COM of mass positions to refined vaules and store COMP
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2
            COMP = [np.round((XCOM)*u.kpc),
                    np.round((YCOM)*u.kpc),
                    np.round((ZCOM)*u.kpc)]

        return COMP

#---Input:  vx,vy,vz and mass
#---Output: returns an array with the x,y,z coords of the COM position

    def COM_V(self, COMX,COMY,COMZ):
        #COMP = self.COM_P(delta)
        #---change to COM frame
        VXNEW = self.x[:]*u.kpc - COMX
        VYNEW = self.y[:]*u.kpc - COMY
        VZNEW = self.z[:]*u.kpc - COMZ
        VRNEW = np.sqrt(VXNEW**2 + VYNEW**2 + VZNEW**2)

        #---store velocities within 8kpc from the COM position
        RVMAX  = 15.0*u.kpc
        vindex = np.where(VRNEW < RVMAX)
        vx2    = self.vx[vindex]
        vy2    = self.vy[vindex]
        vz2    = self.vz[vindex]
        mn2    = self.m[vindex]
        VXCOM, VYCOM, VZCOM = self.COMdefine(vx2,vy2,vz2, mn2)

        #---store and round the COM velocity
        COMV = [np.around((VXCOM)*u.km/u.s),
                np.around((VYCOM)*u.km/u.s),
                np.around((VZCOM)*u.km/u.s)]
        
        return COMV 

#---notes
"""
modified for future assignments
COM_P now takes VolDec as an imput 
      replace RMAX/2 with RMAX/VolDec
"""


#---COM of mass object for MW, M31, M33
MWCOM  = CenterOfMass("MW_000.txt", 2)
M31COM = CenterOfMass("M31_000.txt", 2)
M33COM = CenterOfMass("M33_000.txt", 2)

#---radius and velocity

#--without changing the COMV
#MW_pos = MWCOM.COM_P(1,2)
#print"COM Disk position for MW:",  MW_pos[0], MW_pos[1], MW_pos[2]
#MW_vel = MWCOM.COM_V(1,2)
#print"COM Disk velocity for MW:",  MW_vel[0], MW_vel[1], MW_vel[2]







"""
########################################################
print"Problem 1:"

#---position of MW, M31, M33
MW_pos = MWCOM.COM_P(1,2)
print"COM Disk position for MW:",  MW_pos[0], MW_pos[1], MW_pos[2]
M31_pos = M31COM.COM_P(1,2)
print"COM Disk position for M31:", M31_pos[0], M31_pos[1], M31_pos[2]
M33_pos = M33COM.COM_P(1,2)
print"COM Disk position for M33:", M33_pos[0], M33_pos[1], M33_pos[2]

#---velocity of MW, M31, M33
MW_vel = MWCOM.COM_V(MW_pos[0],MW_pos[1],MW_pos[2])
print"COM Disk velocity for MW:",  MW_vel[0], MW_vel[1], MW_vel[2]
M31_vel = M31COM.COM_V(M31_pos[0],M31_pos[1],M31_pos[2])
print"COM Disk velocity for M31:", M31_vel[0], M31_vel[1], M31_vel[2]
M33_vel = M33COM.COM_V(M33_pos[0],M33_pos[1],M33_pos[2])
print"COM Disk velocity for M33:", M33_vel[0], M33_vel[1], M33_vel[2]
print""


print"Problem 2 & 3:"
#---determine the separation between the MW and M31 
MW_M31 = np.sqrt((M31_pos[0]-MW_pos[0])**2 + (M31_pos[1]-MW_pos[1])**2 + (M31_pos[2]-MW_pos[2])**2)
print("Separation between the MW and M31 =", np.round(MW_M31))

#---determine the relative velocity between the MW and M31
vMW_M31 = np.sqrt((M31_vel[0]-MW_vel[0])**2 + (M31_vel[1]-MW_vel[1])**2 + (M31_vel[2]-MW_vel[2])**2)
print("Relative Velocity between the MW and M31 =", np.round(vMW_M31))

#---determine the separation between M33 and M31
M33_M31 = np.sqrt((M33_pos[0]-M31_pos[0])**2 + (M33_pos[1]-M33_pos[1])**2 + (M33_pos[2]-M31_pos[2])**2)
print("Separation between M33 and M31 =", np.round(M33_M31))

#---determine the relative velocity between M33 and M31
vM33_M31 = np.sqrt((M33_vel[0]-M31_vel[0])**2 + (M33_vel[1]-M31_vel[1])**2 + (M33_vel[2]-M31_vel[2])**2)
print("Relative Velocity between M33 and M31 = ", np.round(vM33_M31))
print""


print"Problem 4:"
print"The iterative process is important because in order to understand the outcome of the merger, one needs to fully know the internal structure of the two objects. You need to know how the MW and M31 behaves in order to make accuate simulations of the merger."

#---pos and vel relative
print"x,y,z pos of m33 realtive to m31"
M33_M31x = (M33_pos[0]-M31_pos[0])
M33_M31y = (M33_pos[1]-M31_pos[1])  
M33_M31z = (M33_pos[2]-M31_pos[2])  
print("x",M33_M31x)
print("y",M33_M31y)
print("z",M33_M31z)

print"vx,vy,vz vel of m33 realtive to m31"
M33_M31vx = (M33_vel[0]-M31_vel[0])
M33_M31vy = (M33_vel[1]-M31_vel[1])  
M33_M31vz = (M33_vel[2]-M31_vel[2])  
print("vx",M33_M31vx)
print("vy",M33_M31vy)
print("vz",M33_M31vz)
"""
