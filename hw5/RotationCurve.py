
"""
Frida Jauregui
Homework 5 started Feb. 8, 2018

Determine the mass distribution and use that to determine each galaxy's roation curve.
"""

import numpy as np
import astropy.units as u
from ReadFile import Read
from CenterofMass import CenterOfMass
import matplotlib.pyplot as plt
from astropy.constants import G
import matplotlib



#---Class
"""
    Follow the structure of the previous assignment
    Input: galaxy: a string with Galaxy name, 
           snap:   snapshot number
"""
class MassProfile:
    def __init__(self, galaxy, Snap):
        #---add string of the filenumber to the vaule "000"
        ilbl = '000' + str(Snap)
        
        #---removing all but last 3 digits
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + '.txt'
        #so if I want the file "M33_000.txt" I would input: MassProfile("M33",00)

        #---read in the data for the x,y,z pos and mass
        self.time, self.total, self.data = Read(self.filename)
        self.m = self.data['m']
        self.x = self.data['x']
        self.y = self.data['y']
        self.z = self.data['z']

        #---store the name of the galaxy as a global property
        self.gname = galaxy
        
#---Functions
        """
        Compute the mass enclosed within a given radius of the COM pos for a specified galaxy
        and a specified galaxy and a specified component of that galaxy
        Input:  particle type
                an array with radii
        Output: returns an array of masses w/units
        """
    def MassEnclosed(self, ptype, R):
        #---determine COM position by creating a com object and calling COM_P
        COM    = CenterOfMass(self.filename, 2)
        GCOM   = COM.COM_P(1.0, 2.0)
        Menc   = np.zeros(np.size(R))
        index  = np.where(self.data['type'] == ptype)
        #---find the COM of the positions
        mG = self.m[index]
        
        xG = self.x[index] - float(GCOM[0]/u.kpc)
        yG = self.y[index] - float(GCOM[1]/u.kpc)
        zG = self.z[index] - float(GCOM[2]/u.kpc)
        
        rG = np.sqrt(xG**2 + yG**2 + zG**2)
        #---loop over the radius to define particles that are enclosed within the radius given
        for i in range(np.size(R)):
            #---store masses within radii
            indexR  = np.where(rG < R[i])
            Menc[i] = np.sum(mG[indexR])
        
        return Menc*1e10*u.Msun
    
        """
        Create a caveat for M33
        Input:  array of radii
        Output: array pf total mass(B+D+H) at each radius of the input array
        """
    def MassEnclosedTotal(self, R):
        #---find mass enclosed for each particle type)
        Menc = self.MassEnclosed(1,R) + self.MassEnclosed(2,R) + self.MassEnclosed(3,R)

        #---for M33 since it contains no bulge
        if (self.gname == 'M33'):
            Menc = self.MassEnclosed(1,R) + self.MassEnclosed(2,R)

        return Menc
    
    """
    Compute the mass enclosed within a given radius
    Input:  radius, scale factor, halo mass
    OutPut: halo mass in Msun units
    """
    def HernquistMass(self, R, a, Mhalo):
        HR = (Mhalo*R**2)/(a+R)**2
        return HR*u.Msun

    """
    Compute using MassEnclosed at each radius
    Input:  ptype, array with radii
    Output: circular speed in units of km/s rounded
    """
    def CircularVelocity(self, ptype, R):
        #---constants
        Gcos = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        Menc = self.MassEnclosed(ptype,R)
        #---rotaion curves V = (G*M(<R)/R)^0.5
        Vcirc = np.round(np.sqrt(G*Menc/R/u.kpc),2)
        
        return Vcirc


    """
    Total CV is NOT just the CV of each individual galaxy component summed together
    Input:  array of radii
    Output: array of CV units of km/s of all galaxy compenents 
            at each radius of the input array
    """
    def CircularVelocityTotal(self, R):
        Menc = self.MassEnclosedTotal(R)
        Vcirc = np.round(np.sqrt(G*Menc/R/u.kpc),2)

        return  Vcirc 

    """
    Compute the circular speed using HR mass profile, can call it here
    Input:  radius, scale factor, halo mass
    Output: circular speed in units of km/s rounded
    """
    def HernquistVCirc(self,R, a, Mhalo):
        Menc = self.HernquistMass(R,a,Mhalo)
        Vcirc = np.round(np.sqrt(G*Menc/R/u.kpc),2)
        #---returns array Vcirc
        
        return Vcirc

#---PLOTS
#---mass profiles
MWP  = MassProfile('MW', 0)
M31P = MassProfile('M31',0)
M33P = MassProfile('M33',0)

#---radius array of 30 kpc
Rarry = np.arange(0.1, 30, 0.5)


########################################################################################
########################################################################################
#---MW
#---plot using inclass lab as template
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
#---halo, disk, bulge, and total of MW
MWH  = MWP.MassEnclosed(1,Rarry)
MWD  = MWP.MassEnclosed(2,Rarry)
MWB  = MWP.MassEnclosed(3,Rarry)
MWT  = MWP.MassEnclosedTotal(Rarry)

#---HR profile want to give total Halo Mass NOT massenc Galaxy Halo Mass
#---using vaules from previous assignments
MWHR  = MWP.HernquistMass(Rarry,  62, 1.975e12)

#---plotting halo, disk, bulge, total, and HR for MW 
plt.semilogy(Rarry, MWH,  color='blue',  label='Halo Mass')
plt.semilogy(Rarry, MWD,  color='red',   label='Disk Mass')
plt.semilogy(Rarry, MWB,  color='green', label='Bulge Mass')
plt.semilogy(Rarry, MWT,  color='yellow',label='Total Mass')
plt.semilogy(Rarry, MWHR, color='black', linestyle="--",label='Hernquist for MW a=62')

#---titles~ things~
plt.xlabel('From Center of Mass [kpc]')
plt.ylabel('Mass Enclosed (Msun)')
plt.ylim(1e9,1e12)
legend = ax.legend(loc='lower right')
plt.figtext(0.15, 0.83, 'Milky Way Mass Profile', fontsize=22)
########################################################################################
#---M31
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

M31H = M31P.MassEnclosed(1,Rarry)
M31D = M31P.MassEnclosed(2,Rarry)
M31B = M31P.MassEnclosed(3,Rarry)
M31T = M31P.MassEnclosedTotal(Rarry)

M31HR = M31P.HernquistMass(Rarry, 62, 1.921e12)

plt.semilogy(Rarry, M31H,  color='blue',  label='Halo Mass')
plt.semilogy(Rarry, M31D,  color='red',   label='Disk Mass')
plt.semilogy(Rarry, M31B,  color='green', label='Bulge Mass')
plt.semilogy(Rarry, M31T,  color='yellow',label='Total Mass')
plt.semilogy(Rarry, M31HR, color='black', linestyle="--",label='Hernquist for M31 a=62')

plt.xlabel('From Center of Mass [kpc]')
plt.ylabel('Mass Enclosed (Msun)')
plt.ylim(1e9,1e12)
legend = ax.legend(loc='lower right')
plt.figtext(0.15, 0.83, 'M31 Mass Profile', fontsize=22)
########################################################################################
#---M33
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

M33H = M33P.MassEnclosed(1,Rarry)
M33D = M33P.MassEnclosed(2,Rarry)
M33T = M33P.MassEnclosedTotal(Rarry)

M33HR = M33P.HernquistMass(Rarry, 25, 0.187e12)

plt.semilogy(Rarry, M33H,  color='blue',  label='Halo Mass')
plt.semilogy(Rarry, M33D,  color='red',   label='Disk Mass')
plt.semilogy(Rarry, M33T,  color='yellow',label='Total Mass')
plt.semilogy(Rarry, M33HR, color='black', linestyle="--",label='Hernquist for M33 a=25')

plt.xlabel('From Center of Mass [kpc]')
plt.ylabel('Mass Enclosed (Msun)')
plt.ylim(1e9,1e12)
legend = ax.legend(loc='lower right')
plt.figtext(0.15, 0.83, 'M33 Mass Profile', fontsize=22)
########################################################################################
########################################################################################
#---rotaion curve MW
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

#---assigning vel for halo mass (disk =2 and bulge =3) and total vel
MWHROT = MWP.CircularVelocity(1,Rarry)
MWDROT = MWP.CircularVelocity(2,Rarry)
MWBROT = MWP.CircularVelocity(3,Rarry)
MWROTT = MWP.CircularVelocityTotal(Rarry)

#---hernquist Vcirc Profile
MWHRCV = MWP.HernquistVCirc(Rarry, 62, 1.975e12)

#---plot halo, disk, bulge, total, CVHR and totCV
plt.plot(Rarry, MWHROT, color='purple', linestyle="--", label='Halo')
plt.plot(Rarry, MWDROT, color='purple', linestyle=":", label='Disk')
plt.plot(Rarry, MWBROT, color='purple', linestyle="-", label='Bulge')

plt.plot(Rarry, MWROTT, color='blue', linestyle="--", label='Total')
plt.plot(Rarry, MWHRCV, color='blue', linestyle="-", label='Hernquist')

plt.xlabel('Radius (kpc)', fontsize=22)
plt.ylabel('Circular Velocity (km/s)', fontsize=22)
#plt.ylim(0,3)
legend = ax.legend(loc='upper right')
plt.figtext(0.15, 0.83, 'Milky Way Circular Velocity', fontsize=22)
########################################################################################
#---rotation curve M31
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

M31HROT = M31P.CircularVelocity(1,Rarry)
M31DROT = M31P.CircularVelocity(2,Rarry)
M31BROT = M31P.CircularVelocity(3,Rarry)
M31ROTT = M31P.CircularVelocityTotal(Rarry)

M31HRCV = M31P.HernquistVCirc(Rarry, 62, 1.921e12)

plt.plot(Rarry, M31HROT, color='purple', linestyle="--", label='Halo')
plt.plot(Rarry, M31DROT, color='purple', linestyle=":", label='Disk')
plt.plot(Rarry, M31BROT, color='purple', linestyle="-", label='Bulge')

plt.plot(Rarry, M31ROTT, color='blue', linestyle="--", label='Total')
plt.plot(Rarry, M31HRCV, color='blue', linestyle="-", label='Hernquist')

plt.xlabel('Radius (kpc)', fontsize=22)
plt.ylabel('Circular Velocity (km/s)', fontsize=22)
legend = ax.legend(loc='upper right')
plt.figtext(0.15, 0.83, 'M31 Circular Velocity', fontsize=22)
########################################################################################
#---rotation curve M33
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

M33HROT = M31P.CircularVelocity(1,Rarry)
M33DROT = M31P.CircularVelocity(2,Rarry)
M33ROTT = M31P.CircularVelocityTotal(Rarry)

M33HRCV = M31P.HernquistVCirc(Rarry, 25, 0.187e12)

plt.plot(Rarry, M33HROT, color='purple', linestyle="--", label='Halo')
plt.plot(Rarry, M33DROT, color='purple', linestyle=":", label='Disk')

plt.plot(Rarry, M33ROTT, color='blue', linestyle="--", label='Total')
plt.plot(Rarry, M33HRCV, color='blue', linestyle="-", label='Hernquist')

plt.xlabel('Radius (kpc)', fontsize=22)
plt.ylabel('Circular Velocity (km/s)', fontsize=22)
legend = ax.legend(loc='upper right')
plt.figtext(0.15, 0.83, 'M33 Circular Velocity', fontsize=22)

plt.show()
   
        




