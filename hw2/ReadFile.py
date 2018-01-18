"""
Frida Jauregui
ASTR400b HW2 Due Jan 23
"""

import numpy as np
import astropy.units as u

def Read(MW_000):
    file = open(MW_000, 'r')         #open text file
    line1 = file.readline()          #reading first line and storing the time
    label, value = line1.split()
    time = float(value)*10.0*u.Myr

    line2 = file.readline()
    label, value = line2.split()     #reading first line and storing particles
    particles = float(value)

    file.close()

    data = np.genfromtxt(MW_000, dtype=None, names=True, skip_header=3)
    #storing the remainder of the MW_000 file to be able to call back
    return time, particles, data



    
    

    
    
