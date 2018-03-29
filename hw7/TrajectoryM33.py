"""
Frida Jauregui
Homework 7 started March 5, 2018
predict the trajectory of M33 as it orbits M31.
use the COM file for M31 and M33 from hw6
"""

import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib

from ReadFile import Read
from CenterofMass import CenterOfMass

#--Class
#---create a series of functions to determine the acceleration M33 feels from M31
class M33AnalyticOrbit:
    def __init__(self, filename):
        #--G constant in kpc^3/Gyr^2/Msun
        self.G = 4.498768e-6
        ###self.time, self.total, self.data = Read(filename)
        self.outfile = filename
        #--com pos.&vel. of M33 relative to M31 using disk particles for snapshot 0
	self.x  = -98.0
        self.y  = -120.0
        self.z  = -127.0
        self.vx = -29.0
        self.vy = 174.0
        self.vz = 93.0
        #--use Hernquist scale length from assignment 5 and masses from assignment 3
        self.rd    = 5.0
        self.Mdisk = 0.12e12
        
        self.rbulge = 1.0
        self.Mbulge = 0.019e12
        
        self.rhalo = 62.0
        self.Mhalo = 1.921e12

        
#--Functions
#---define functions that will compute the gravitational acceleration from 3 comp. of M31
    def HenquistAccel(self, M, ra, x,y,z, ddv):
#Input: M is total halo or bulge mass
#------ ra is the scale length
#------ x,y,z coord. of the gravitaional accel.
#------ av is a dummy varrible that indicates which comp. of the accel.
#Output:returns accel. in the direction of the input of the dummy varible
   
        #--grav. accel. induced by a hernquist profile
        r  = np.sqrt(x**2 + y**2 + z**2)
        ###ax = -G*M*x/((r(ra + r))**2)
        ###ay = -G*M*y/((r(ra + r))**2)
        ###az = -G*M*z/((r(ra + r))**2)
        dom = r*(ra+r)**2
        #--set a dummy variable that indicates if x, y, or z
        if ddv =='x':
            av = x
        if ddv =='y':
            av = y
        if ddv =='z':
            av = z
        ah  = -self.G*M*av/dom
        return ah

    
    def MiyamotoNagaiAccel(self, M, rd, x,y,z, ddv):
#Input: M is total disk mass
#------ rd is the scale length, disk
#------ x,y,z coord. of the gravitaional accel.
#------ av is a dummy varrible that indicates which comp. of the accel.
#Output:returns accel. in the direction of the input of the dummy varible
    
        #--define zd   
        zd = self.rd/5.0
        #--varibles from Miyamoto-Nagai 1975 profile
        R  = np.sqrt(x**2 + y**2)
        B  = self.rd + np.sqrt(z**2 + zd**2)

        adx = -self.G*M*x/((R**2 + B**2)**1.5)
        ady = -self.G*M*y/((R**2 + B**2)**1.5)
        #--accel. in the z direction is different
        RBdom = (R**2 + B**2)**(1.5)
        Zdom  = np.sqrt(z**2 + zd**2)
        adz = (-self.G*M*B*z)/(RBdom*Zdom)
        #--dummy variable
        if ddv =='x':
            am = adx
        if ddv =='y':
            am = ady
        if ddv =='z':
            am = adz
        return am
    

    def M31Accel(self, x,y,z, ddv):
#Input: x,y,z coord. of the gravitaional accel.
#------ ddv is a dummy varrible that indicates which comp. of the accel.
#Output:sums all acceleration terms from each galaxy component (disk,bluge,halo)
        
        haloaccel  = self.HenquistAccel(self.Mhalo, self.rd, x,y,z, ddv)
        bulgeaccel = self.HenquistAccel(self.Mbulge, self.rbulge, x,y,z, ddv)
        diskaccel  = self.MiyamotoNagaiAccel(self.Mdisk, self.rhalo, x,y,z, ddv)
        allaccels  = haloaccel + bulgeaccel + diskaccel
        return allaccels

    
    def LeapFrog(self, dt, x,y,z, vx,vy,vz):
#Input: dt is a time interval for integration
#------ x,y,z is a starting potential for the M33 COM pos.
#------ starting potential for the M33 COM vel.
#Output:array of pos&vel

        #--bulid an integrator
        xn = x + (dt/2)*vx
        yn = y + (dt/2)*vy
        zn = z + (dt/2)*vz

        M31anx = self.M31Accel(xn,yn,zn,'x')
        M31any = self.M31Accel(xn,yn,zn,'y')
        M31anz = self.M31Accel(xn,yn,zn,'z')
        
        vxn = vx +  M31anx*dt
        vyn = vy +  M31any*dt
        vzn = vz +  M31anz*dt

        rx = x + 0.5*(vx + vxn)*dt
        ry = y + 0.5*(vy + vyn)*dt
        rz = z + 0.5*(vz + vzn)*dt
        return rx,ry,rz, vxn,vyn,vzn

    
    def OrbitIntegrator(self, t0, dt, tmax):
#Input: t0 starting time
#------ dt a time interval
#------ tmax final time
#Output:array of of motion

        #--supply starting COM pos&vel of M33 realtive to M31
        x = self.x
        y = self.y
        z = self.z

        vx = self.vx
        vy = self.vy
        vz = self.vz

                
        #--initialize the array
        a = int((tmax)/dt)+1
        OrbitM33 = np.zeros((a,7))

        t = t0
        #--while loop over leapfrog
        while t < tmax:
            OrbitM33[int(t/dt),0] = t
            OrbitM33[int(t/dt),1] = x
            OrbitM33[int(t/dt),2] = y
            OrbitM33[int(t/dt),3] = z
            OrbitM33[int(t/dt),4] = vx
            OrbitM33[int(t/dt),5] = vy
            OrbitM33[int(t/dt),6] = vz
            
            #--new pos&vel
            x,y,z,vx,vy,vz = self.LeapFrog(dt, x, y, z, vx, vy, vz)
            #--updating time
            t = t+dt

        np.savetxt('TrajectoryM33.txt', OrbitM33, header='t   x    y    z   vx   vy   vz',
                   comments='# ', fmt=['%.2f', '%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])
        
#---plots
M33ana  = M33AnalyticOrbit('TrajectoryM33.txt')
intM33  = M33ana.OrbitIntegrator(0,0.1,10)

#--reading m31 and m33 orbit files from hw6
M31 = np.genfromtxt('M31_lifetime.txt',dtype=None,names=True)
M33 = np.genfromtxt('M33_lifetime.txt',dtype=None,names=True)
M33_trajectory = np.genfromtxt('TrajectoryM33.txt',dtype=None,names=True)

#--plot POS&VEL of M33 relative to M31
time = M33_trajectory['t'] 
POS = np.sqrt(M33_trajectory['x']**2 + M33_trajectory['y']**2 + M33_trajectory['z']**2)
M31time = M31['t']
xsep = (M31['x'] - M33['x'])
ysep = (M31['y'] - M33['y'])
zsep = (M31['z'] - M33['z'])
mag  = np.sqrt(xsep**2 + ysep**2 + zsep**2)

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.plot(time,POS,label='Analytic')
plt.plot(M31time,mag,label='Simulation')
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Separation (kpc)', fontsize=22)
plt.title('M31 and M33 Separation', fontsize=22)
plt.legend()


VEL = np.sqrt(M33_trajectory['vx']**2 + M33_trajectory['vy']**2 + M33_trajectory['vz']**2)
vxsep = (M31['vx'] - M33['vx'])
vysep = (M31['vy'] - M33['vy'])
vzsep = (M31['vz'] - M33['vz'])
mag2  = np.sqrt(xsep**2 + ysep**2 + zsep**2)

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.plot(time,VEL,label='Analytic')
plt.plot(M31time,mag2,label='Simulation')
plt.xlabel('Time (Gyr)', fontsize=22)
plt.ylabel('Separation (km/s)', fontsize=22)
plt.title('M31 and M33 Velocity', fontsize=22)
plt.legend()
plt.show()

#---problems
print("question 1:")
print("So far, I have zeros in my data set and my graphs are all wrong")
print("")
print("question 2:")
print("What's missing the are the tidal stripping effects M33 goes through when it passes through M31")
print("")
print("question 3:")
print("To add the effects MW might have a similar approach might be needed, like with M33ana, consider MW as point mass.")



