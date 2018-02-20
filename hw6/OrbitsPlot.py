"""
Frida Jauregui
Homework 6 started Feb. 16, 2018

store the seperation and relative vel of MW, M31 & M33 over 
the entire simulation and plot the corresponding orbits

compute orbits up to snapshot 800(~12Gyr) of intervals of n=5

"""

import numpy as np
import astropy.units as u
from ReadFile import Read
from CenterofMass import CenterOfMass
from Orbits import OrbitCOM
import matplotlib.pyplot as plt

#--compute the time, comp and comv for each galaxy from snap 0 to 800
#--generate text files
MW  = np.genfromtxt('MW_lifetime.txt', dtype = float, names=True)
M31 = np.genfromtxt('M31_lifetime.txt', dtype = float, names=True)
M33 = np.genfromtxt('M33_lifetime.txt', dtype = float, names=True)

#--MW values
MWtime = MW['t']
MWx    = MW['x']
MWy    = MW['y']
MWz    = MW['z']
#print(MWz)

MWvx = MW['vx']
MWvy = MW['vy']
MWvz = MW['vz']

#--M31values
M31time = M31['t']
M31x    = M31['x']
M31y    = M31['y']
M31z    = M31['z']

M31vx = M31['vx']
M31vy = M31['vy']
M31vz = M31['vz']

#--M33 values
M33time = M33['t']
M33x    = M33['x']
M33y    = M33['y']
M33z    = M33['z']

M33vx = M33['vx']
M33vy = M33['vy']
M33vz = M33['vz']

#--COM seperation
xsep = M31x - MWx
ysep = M31y - MWy
zsep = M31z - MWz
mag  = np.sqrt(xsep**2 + ysep**2 + zsep**2)

xsep2 = M31x - M33x
ysep2 = M31y - M33y
zsep2 = M31z - M33z
mag2  = np.sqrt(xsep2**2 + ysep2**2 + zsep2**2)

#--Vel sep
vxsep = M31vx - MWvx
vysep = M31vy - MWvy
vzsep = M31vz - MWvz
vmag  = np.sqrt(vxsep**2 + vysep**2 + vzsep**2)

vxsep2 = M31vx - M33vx
vysep2 = M31vy - M33vy
vzsep2 = M31vz - M33vz
vmag2  = np.sqrt(vxsep2**2 + vysep2**2 + vzsep2**2)

###POS
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('Seperation Between M31 and MW')
plt.xlabel('Time [Gyr]')
plt.ylabel('Magnitude of Seperation [kpc]')
plt.plot(MW['t'], mag,  color='red')

ax.set_rasterized(True)
#plt.savefig('M31_MWpossep.png', rasterized=True, dpi=350)


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('Seperation Between M31 and M33')
plt.xlabel('Time [Gyr]')
plt.ylabel('Magnitude of Seperation [kpc]')
plt.plot(M33['t'], mag2,  color='blue')

ax.set_rasterized(True)
#plt.savefig('M31_M33possep.png', rasterized=True, dpi=350)



###VEL
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('Velocity Seperation Between M31 and MW')
plt.xlabel('Time [Gyr]')
plt.ylabel('Magnitude of Vel. Seperation [km/s]')
plt.plot(MW['t'], vmag,  color='green')

ax.set_rasterized(True)
#plt.savefig('M31_MWvelsep.png', rasterized=True, dpi=350)


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.title('Velocity Seperation Between M31 and M33')
plt.xlabel('Time [Gyr]')
plt.ylabel('Magnitude of Vel. Seperation [km/s]')
plt.plot(M33['t'], vmag2,  color='black')

ax.set_rasterized(True)
#plt.savefig('M31_M33velsep.png', rasterized=True, dpi=350)


plt.show()

###QUESTIONS
print("Question 1:")
print("As the magnitude of seperation reaches zero it will be a close encounter, so before it stays at zero MW and M31 will experience two close encounters.")
print("")

print("Question 2:")
print("Viewing the velocity seperation and the position seperation of M31 and MW I can see that as the two galaxies experience a close encounter, their velocities will increase, which makes sense considering close encounters with objects will cause the objects to zip past each other at a great speed.")
print("")

print("Question 3:")
print("I would say around ~6.3Gyr would be the time for M31 and M31 to merge. At around this time the velocity seperation of M31 and M33 is at its lowest while its position seperation is slowly decreasing.")
print("")


