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
import pandas as pd

from ReadFile import Read
from CenterofMass import CenterOfMass

G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

#--Class
#---create a series of functions to determine the acceleration M33 feels from M31
class M33AnalyticOrbit:
    def __init__(self, filename):
        self.time, self.total, self.data = Read(filename)
	self.x  = -98.0*u.kpc
        self.y  = -120.0*u.kpc
        self.z  = -127.0*u.kpc
        self.vx = -29.0*(u.km/u.s)
        self.vy = -174.0*(u.km/u.s)
        self.vz = 93.0*(u.km/u.s)
        #--com pos.&vel. of M33 relative to M31 using disk particles for snapshot 0
        self.pM33_M31 = 160*u.kpc
        self.vM33_M31 = 199.0*(u.km/u.s)
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

        #adx = -G*M*x/((R**2 + B**2)**1.5)
        #ady = -G*M*y/((R**2 + B**2)**1.5)
        #--accel. in the z direction is different
        adz = -G*M*B*z/((((R**2 + B**2)**1.5)*(z**2 + zd**2)**0.5))
        #--dummy variable
        am = (-G*M*ddv/(R**2 + B**2)**1.5)
        if (av == z):
            a2 = adz
        return am 



    def M31Accel(x,y,z, ddv):
#Input: x,y,z coord. of the gravitaional accel.
#------ av is a dummy varrible that indicates which comp. of the accel.
#Output:sums all acceleration terms from each galaxy component (disk,bluge,halo)

        haloaccel = self.HenquistAccel(1.921*10**12, 62, x,y,z, ddv)
        bulgaccel = self.HenquistAccel(0.019*10**12, 62, x,y,z, ddv)
        diskaccel = self.MiyamotoNagaiAccel(0.12*10**12, 62, x,y,z, ddv)
        allaccels = haloaccel + bulgeaccel + diskaccel
        return allacccels

        
    N = 1000
    t = np.linspace(0,100,N)
    dt = t[1] - t[0]
    def LeapFrog(dt, x,y,z, vx,vy,vz):
#Input: dt is a time interval for integration
#------ x,y,z is a starting potential for the M33 COM pos.
#------ starting potential for the M33 COM vel.
#Output:equation of motion

    #--arrays filled with zeros
        x = np.zeros(N)
        y = np.zeros(N)
        z = np.zeros(N)
        vx = np.zeros(N)
        vy = np.zeros(N)
        vz = np.zeros(N)
    
        for i in range(1,N):
            an = self.M31Accel(x,y,z,ddv)
            x[i+0.5] = x[i] + (dt/2)*vx[i]
            y[i+0.5] = x[i] + (dt/2)*vy[i]
            z[i+0.5] = x[i] + (dt/2)*vz[i]
        
            vx[i+1] = vx[i] + an[i+0.5]*dt
            vy[i+1] = vy[i] + an[i+0.5]*dt
            vz[i+1] = vz[i] + an[i+0.5]*dt

            rx[i+1] = rx[i] + 0.5*(vx[i] + vx[i+1])*dt
            ry[i+1] = ry[i] + 0.5*(vy[i] + vy[i+1])*dt
            rz[i+1] = rz[i] + 0.5*(vz[i] + vz[i+1])*dt
            
        mot = np.sqrt(rx**2 + ry**2 + rz**2)
        return mot

    def OrbitIntegrator(self, t0, dt, tmax):
        x = self.x
        y = self.y
        z = self.z

        vx = self.vx
        vy = self.vy
        vz = self.vz
