"""
Frida Jauregui
Homework 7 started March 5, 2018
predict the trajectory of M33 as it orbits M31.
use the COM file for M31 and M33 from hw6
"""

import numpy as np
import astropy.units as u
from astropy.constants import G
import matplotlib.pyplot as plt
import matplotlib

from ReadFile import Read
from CenterofMass import CenterOfMass

G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

#--Class
#---create a series of functions to determine the acceleration M33 feels from M31
class M33AnalyticOrbit:
    def __init__(self, filename):
        self.time, self.total, self.data = Read(filename)
        #--com pos.&vel. of M33 relative to M31 using disk particles for snapshot 0
	self.x  = -98.0*u.kpc
        self.y  = -120.0*u.kpc
        self.z  = -127.0*u.kpc
        self.vx = -29.0*(u.km/u.s)
        self.vy = -174.0*(u.km/u.s)
        self.vz = 93.0*(u.km/u.s)
        #--use Hernquist scale length from assignment 5 and masses from assignment 3
        self.rd    = 5*u.kpc
        self.Mdisk = (0.12*10**12)*u.Msun
        
        self.rbulge = 1.0*u.kpc
        self.Mbulge = (0.019*10**12)*u.Msun
        
        self.rhalo = 62*u.kpc
        self.Mhalo = (1.921*10**12)*u.Msun

        
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
        ax = -G*M*x/((r(ra + r))**2)
        ay = -G*M*y/((r(ra + r))**2)
        az = -G*M*z/((r(ra + r))**2)
        #--set a dummy variable that indicates if x, y, or z
        ddv = '%s_'%(x,y,z)
        ah  = -G*M*ddv/((r(ra + r))**2)
        return ah

    
    def MiyamotoNagaiAccel(self, M, rd, x,y,z, ddv):
#Input: M is total disk mass
#------ rd is the scale length, disk
#------ x,y,z coord. of the gravitaional accel.
#------ av is a dummy varrible that indicates which comp. of the accel.
#Output:returns accel. in the direction of the input of the dummy varible
    
        #--define zd    
        zd = (self.rd)/5.0
        #--varibles from Miyamoto-Nagai 1975 profile
        R  = np.sqrt(x**2 + y**2)
        B  = rd + np.sqrt(z**2 + zd**2)

        adx = -G*M*x/((R**2 + B**2)**1.5)
        ady = -G*M*y/((R**2 + B**2)**1.5)
        #--accel. in the z direction is different
        adz = -G*M*B*z/((((R**2 + B**2)**1.5)*(z**2 + zd**2)**0.5))
        #--dummy variable
        ddv = '%s_'%(x,y,z)
        am = (-G*M*ddv/(R**2 + B**2)**1.5)
        if (ddv == z):
            am = adz
        return am
    

    def M31Accel(x,y,z, ddv):
#Input: x,y,z coord. of the gravitaional accel.
#------ ddv is a dummy varrible that indicates which comp. of the accel.
#Output:sums all acceleration terms from each galaxy component (disk,bluge,halo)
        
        haloaccel = self.HenquistAccel(1.921*10**12, 62, x,y,z, ddv)
        bulgaccel = self.HenquistAccel(0.019*10**12, 62, x,y,z, ddv)
        diskaccel = self.MiyamotoNagaiAccel(0.12*10**12, 62, x,y,z, ddv)
        allaccels = haloaccel + bulgeaccel + diskaccel
        return allacccels

    
    def LeapFrog(self, dt, x,y,z, vx,vy,vz):
#Input: dt is a time interval for integration
#------ x,y,z is a starting potential for the M33 COM pos.
#------ starting potential for the M33 COM vel.
#Output:equation of motion

        #--bulid an integrator
        ddv = '%s_'%(x,y,z)
        an = self.M31Accel(x,y,z,ddv)
        x = x + (dt/2)*vx
        y = y + (dt/2)*vy
        z = z + (dt/2)*vz
        
        vx = vx + an*dt
        vy = vy + an*dt
        vz = vz + an*dt

        rx = x + 0.5*(vx + vx)*dt
        ry = y + 0.5*(vy + vy)*dt
        rz = z + 0.5*(vz + vz)*dt
        return vx,vy,vz, rx,ry,rz

    
    def OrbitIntegrator(self, t0, dt, tmax):
#Input: t0 starting time
#------ dt a time interval
#------ tmax final time
#Output:array of of motion

        fileout = 'TrajectoryM33.txt'
        #--supply starting COM pos&vel of M33 realtive to M31
        x = self.x
        y = self.y
        z = self.z

        vx = self.vx
        vy = self.vy
        vz = self.vz

        #--initialize the array
        a = int(tmax/dt)+1
        OrbitM33 = np.zeros((a,7))  
        t  = t0
        #--while loop over leapfrog
        while(t < tmax):
            vx,vy,vz,x,y,z = self.LeapFrog(dt, x,y,z, vx,vy,vz)

            OrbitM33[int(tmax/dt),0] = t/u.Myr/1000
            OrbitM33[int(tmax/dt),1] = x
            OrbitM33[int(tmax/dt),2] = y
            OrbitM33[int(tmax/dt),3] = z
            OrbitM33[int(tmax/dt),4] = vx
            OrbitM33[int(tmax/dt),5] = vy
            OrbitM33[int(tmax/dt),6] = vz

        np.savetxt(fileout, OrbitM33, header='t   x    y    z   vx   vy   vz', comments='# ',
                fmt=['%.2f', '%.2f','%.2f','%.2f','%.2f','%.2f','%.2f'])

print("trying out code")
M33 =  M33AnalyticOrbit('M33_000.txt')
om33 = M33.OrbitIntegrator(0,0.5,10)
