"""
Frida Jauregui
Started Jan 16
ASTR400b HW2 Due Jan 23
"""

import numpy as np
import astropy.units as u

#---Data
data = np.genfromtxt("MW_000.txt",dtype=None,names=True,
                     skip_header=3)

#---Want Disk particle
indextype = np.where(data['type']>1)
newdata = data[:][indextype]
"""
#---check
print"newdata:"
print(newdata)
"""

#---Mass
Mindex = np.where(newdata['m'][:100])             #the mass of the 100th particle
mnew = newdata['m'][Mindex]*1.0e10*u.solMass      #assigning mass from data in the text

#---Position            
Xindex = np.where(newdata['x'][:100])             #the x position of the 100th particle
dxnew = newdata['x'][Xindex]*1.0*u.kpc            #assigning pos and quantity
                      
Yindex = np.where(newdata['y'][:100])             #the y position of the 100th particle
dynew = newdata['y'][Yindex]*1.0*u.kpc
    
Zindex = np.where(newdata['z'][:100])             #the z position of the 100th particle
dznew = newdata['z'][Zindex]*1.0*u.kpc
    
#---Velocity
Vxindex = np.where(newdata['vx'][:100])           #the x velocity of the 100th particle
vxnew = newdata['vx'][Vxindex]*1.0*u.kilometer/u.second
        
Vyindex = np.where(newdata['vy'][:100])           #the y velocity of the 100th particle
vynew = newdata['vy'][Vyindex]*1.0*u.kilometer/u.second         

Vzindex = np.where(newdata['vz'][:100])           #the z velocity of the 100th particle
vznew = newdata['vz'][Vzindex]*1.0*u.kilometer/u.second

#---3D distance
distance = np.sqrt(dxnew[-1]**2 + dynew[-1]**2 + dznew[-1]**2)
print"3D-distance:"
print(distance)

#---3D
velocity = np.sqrt(vxnew[-1]**2 + vynew[-1]**2 + vznew[-1]**2)
print"3D-velocity:"
print(velocity)

#---convert distance using astropy
lyr = u.lyr
xconv = dxnew.to(lyr)
yconv = dynew.to(lyr)
zconv = dznew.to(lyr)
D = distance.to(lyr) 

print"Mass, X Position, Y Position, Z Position, X Velocity, Y Velocity, Z Velocity"
print(np.around(mnew[-1],3), np.around(xconv[-1],3),
      np.around(yconv[-1],3), np.around(zconv[-1],3),
      np.around(vxnew[-1],3), np.around(vynew[-1],3),np.around(vznew[-1],3))
print"3D-distance Converted to Lightyears:"
print(np.around(D,3))

#---notes
"""
np.where(data['dtype'][:100]) gives me 100 dtype data points, I can print the
last line of the array by array[-1]

keeping this as a function I couldnt print anything back, it just compiled so I just wrote
this assignment as lines of code
"""

    
    

