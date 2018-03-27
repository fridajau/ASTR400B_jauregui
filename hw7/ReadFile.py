"""
Frida Jauregui
Started Jan 16
ASTR400b HW2 Due Jan 23
"""
import numpy as np
import astropy.units as u

#---open and read txt file
def Read(filename):
    file = open(filename,'r')              #open the file and read
    line1 = file.readline()                #read line 1
    label, value = line1.split()           
    time = float(value)*u.Myr              #store line 1 units

    line2 = file.readline()                
    label, value = line2.split()           #assign line 2 a label and value
    total  = float(value)

    file.close()                           #close file
    
    #---store the rest of the file starting from line 4, print/check data
    data = np.genfromtxt(filename,dtype=None,names=True,
                     skip_header=3)

    #---Return
    return time, total, data

"""
#---try to extract  data and understand

data = np.genfromtxt("MW_000.txt",dtype=None,names=True,
                     skip_header=3)
index = np.where(data['x'][:100])
xnew = data['x'][index]
print"X-position:"
print(xnew[-1])
"""

#---notes
"""
revised for hw6
took away a factor of 10 in storing the time
"""
