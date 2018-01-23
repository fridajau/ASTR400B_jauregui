"""
Frida Jauregui
Started Jan 16
ASTR400b HW2 Due Jan 23
"""

import numpy as np
import astropy.units as u

data = np.genfromtxt("MW_000.txt",dtype=None,names=True,
                     skip_header=3)
#---Mass
Mindex = np.where(data['m'][:100])             #the mass of the 100th particle
mnew = data['m'][Mindex]*1.0e10*u.solMass      #assigning mass from data in the text

#---Position            
Xindex = np.where(data['x'][:100])             #the x position of the 100th particle
dxnew = data['x'][Xindex]*1.0*u.kpc            #assigning pos and quantity
                      
Yindex = np.where(data['y'][:100])             #the y position of the 100th particle
dynew = data['y'][Yindex]*1.0*u.kpc
    
Zindex = np.where(data['z'][:100])             #the z position of the 100th particle
dznew = data['z'][Zindex]*1.0*u.kpc
    
#---Velocity
Vxindex = np.where(data['vx'][:100])           #the x velocity of the 100th particle
vxnew = data['vx'][Vxindex]*1.0*u.kilometer/u.second
        
Vyindex = np.where(data['vy'][:100])           #the y velocity of the 100th particle
vynew = data['vy'][Vyindex]*1.0*u.kilometer/u.second         

Vzindex = np.where(data['vz'][:100])           #the z velocity of the 100th particle
vznew = data['vz'][Vzindex]*1.0*u.kilometer/u.second

#---convert distance using astropy
kms = u.kilometer / u.second
lyr = u.lyr / u.yr
kms.to(lyr)


print"Mass, X Position, Y Position, Z Position, X Velocity, Y Velocity, Z Velocity"
print(mnew[-1],dxnew[-1],dynew[-1],dznew[-1],vxnew[-1],vynew[-1],vznew[-1])



#---notes
"""
np.where(data['dtype'][:100]) gives me 100 dtype data points, I can print the
last line of the array by array[-1]

keeping this as a function I couldnt print anything back, it just compiled so I just wrote
this assignment as lines of code
"""

    
    

