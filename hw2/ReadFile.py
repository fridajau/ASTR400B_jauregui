"""
Frida Jauregui
Started Jan 16
ASTR400b HW2 Due Jan 23
"""
import numpy as np
import astropy.units as u


filename = "MW_000.txt"                #assign the file
file = open(filename,'r')              #open the file and read
line1 = file.readline()                #read line 1
label, value = line1.split()           
time = float(value)*10.0*u.Myr         #store line 1 units

line2 = file.readline()                
label, value = line2.split()           #assign line 2 a label and value
particles = float(value)

print"Time, Particles"
print(time, particles)

file.close()                           #close file
    
#store the rest of the file starting from line 4
data = np.genfromtxt("MW_000.txt",dtype=None,names=True,
                     skip_header=3)
#return time, total particles and data to call back
print"Data:"
print(data)

"""
###----check and extract  data

data = np.genfromtxt("MW_000.txt",dtype=None,names=True,
                     skip_header=3)
index = np.where(data['x'][:100])
xnew = data['x'][index]
print"X-position:"
print(xnew[-1])
"""

###---notes
"""
keeping this as a function I couldnt print anything back, it just 
compiled so I just wrote this assignment as lines of code and printed 
the time, total number of particles, and the data as lines of code 
"""
