"""
Frida Jauregui
ASTR400b HW2 Due Jan 23
"""

import numpy as np
import astropy.units as u
from ReadFile import Read 

def ParticleInfo(distance,velocity,mass):
    np.around(vaule,3)
    mass = data['m']*1.0e10*u.M_solarmass          #assigning mass from data in ReadFile
    index = np.where(data['m'][:100])              #the mass of the 100th particle
    mnew = data['m'][index]                        #mass of the particle with the given property
    distance = data['x','y','z']*u.kpc
    iindex = np.where(data['x','y','z'][:100])     #the distance of the 100th particle
    dnew = data['x','y','z'][iindex]
    velocity = data['vx','vy','vz']*u.km/s
    iiindex = np.where(data['vx','vy','vz'][:100]) #the velocity of the 100th particle
    vnew = data['vx','vy','vz'][iiindex]

    return dnew,vnew,mnew
    
    
    
