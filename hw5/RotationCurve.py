
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
        COM    = CenterOfMass(self.filename, ptype)
        COMpos = COM.COM_P(1.0)
        Mass   = np.zeros(np.size(R))
        index  = np.where(self.data['type'] == ptype)
        #---find the COM of the positions
        m2 = self.m[index]
        
        x2 = self.x[index] - float(COMpos[0]/u.kpc)
        y2 = self.y[index] - float(COMpos[1]/u.kpc)
        z2 = self.z[index] - float(COMpos[2]/u.kpc)
        
        R2 = np.sqrt(x2**2 + y2**2 + z2**2)
        #---loop over the radius to define particles that are enclosed within the radius given
        for i in range(np.size(R)):
            #---store masses within radii
            index2  = np.where(R2 < R[i])
            mnew    = m2[index2]
            #---mass enclosed
            massenc = np.sum(mnew)
            #---returns an array
            Mass[i] = (massenc)
        
        return Mass*1e10*u.Msun
    
        """
        Create a caveat for M33
        Input:  array of radii
        Output: array pf total mass(B+D+H) at each radius of the input array
        """
    def MassEnclosedTotal(self, R):
        #---find mass enclosed for each particle type
        DiskMass = self.MassEnclosed(1,R)
        HaloMass = self.MassEnclosed(2,R) 

        #---for M33 since it contains no bulge
        if (self.gname == 'M33'):
            BulgeMass = np.zeros(np.size(R))
        else:
            BulgeMass = self.MassEnclosed(3,R)
            
        BulgeMass= self.MassEnclosed(3,R)
        masstot  = (DiskMass + HaloMass + BulgeMass)
        
        return masstot
    
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
    def CircularVelocity(self, pytpe, R):
        #---constants
        Gcos = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        #---like in mass enclosed follow the loop
        #---rotaion curves V = (G*M(<R)/R)^0.5
        #---loop over the radius
        Vel = np.zeros(np.size(R))
        num = 0
        for i in range(np.size(R)):
            v      = np.sqrt(Gcos*self.MassEnclosed(ptype,R)[num]/radii)
            Vel[i] = v
            num    = num+1
            
        return np.around(Vel,2)*u.km/u.s

    """
    Total CV is NOT just the CV of each individual galaxy component summed together
    Input:  array of radii
    Output: array of CV units of km/s of all galaxy compenents 
            at each radius of the input array
    """
    def CircularVelocityTotal(self, R):
        return R

    """
    Compute the circular speed using HR mass profile, can call it here
    Input:  radius, scale factor, halo mass
    Output: circular speed in units of km/s rounded
    """
    def HernquistVCirc(self,R):
        return R

#---PLOTS

#---plot using inclass lab as template
figure = plt.figure(1)
ax = plt.subplot(111)

#---radius array of 30 kpc
Rarry = np.arange(0.1, 30, 0.5)

#---mass profiles
MWP  = MassProfile('MW', 00)
M31P = MassProfile('M31', 00)
M33P = MassProfile('M33', 00)

#---halo, disk, bulge, and total of MW
MWH  = MWP.MassEnclosed(1,Rarry)
M31H = M31P.MassEnclosed(1,Rarry)
M33H = M33P.MassEnclosed(1,Rarry)
MWD  = MWP.MassEnclosed(2,Rarry)
MWB  = MWP.MassEnclosed(3,Rarry)
MWT  = MWP.MassEnclosedTotal(Rarry)

#print("MWH:", MWH)
#print("MWD:", MWD)
#print("MWB:", MWB)

#---HR profile want to give total Halo Mass NOT massenc Galaxy Halo Mass
#---had help from Drew about the a vaules, since I was running out of time
MWHR  = MWP.HernquistMass(Rarry,  62, MWH)
M31HR = M31P.HernquistMass(Rarry, 62, M31H)
M33HR = M33P.HernquistMass(Rarry, 25, M33H)

#---plotting halo, disk, bulge, total, and HR 
plt.semilogy(Rarry, MWH,  color='blue',  label='Halo Mass')
plt.semilogy(Rarry, MWD,  color='red',   label='Disk Mass')
plt.semilogy(Rarry, MWB,  color='green', label='Bulge Mass')
plt.semilogy(Rarry, MWT,  color='yellow',label='Total Mass')

plt.semilogy(Rarry, MWHR, color='black', linestyle="--",label='Hernquist for MW a=62')
plt.semilogy(Rarry, M31HR, color='black', linestyle=":",label='Hernquist for M31 a=62')
plt.semilogy(Rarry, M33HR, color='black', linestyle="-",label='Hernquist for M33 a=25')

#---titles~ things~
plt.title('Milky Way Mass Profile')
plt.xlabel('From Center of Mass [kpc]')
plt.ylabel('Mass Enclosed (Msun)')
legend = ax.legend(loc='lower right')


#############rotation curve################
#---plotting rotaion curve
figure = plt.figure(2)
a2x = plt.subplot(111)

MWROT  =  MWP.CircularVelocity(1,Rarry)
M31ROT = M31P.CircularVelocity(1,Rarry)
M33ROT = M33P.CircularVelocity(1,Rarry)

plt.semilogy(Rarry, MWROT, color='blue', linestyle="--", label='Halo Mass for MW')
plt.semilogy(Rarry, M31ROT, color='blue', linestyle=":", label='Halo Mass for M31')
plt.semilogy(Rarry, M33ROT, color='blue', linestyle="-", label='Halo Mass for M33')

#---titles~ things~
plt.title('Milky Way Rotation Curve')
plt.xlabel('From Center of Mass [kpc]')
plt.ylabel('Mass Enclosed (Msun)')
legend2 = a2x.legend(loc='lower right')
plt.show()

#---notes
"""
I get an error from plotting the rotation curve
 File "RotationCurve.py", line 121, in CircularVelocity
    v= np.sqrt(Gcos*self.MassEnclosed(ptype,R)[num]/radii)
    NameError: global name 'ptype' is not defined
"""



        
        
        




