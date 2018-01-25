"""
Frida Jauregui
ASTR400b HW3 Due Jan 30th
Started Jan 25th

compute the mass breakdown of the Local Group (snap# 0) using MW, M31, and M31
"""

import numpy as np
import astropy.units as u
from ReadFile import Read

"""
#---notes
input PType:    Halo(1), Disk(2), Bluge(3) 
      filename: "MW00,MW31,MW33"
return total mass
"""
def ComponentMass(PType, filename):
    #---to be read in the file
    time, total, data = Read(filename)

    #---an array to store Ptype
    index = np.where(data['type'] == PType)

    #---Mass
    mnew = data['m'][index]*1e10*u.Msun
    totalmass = np.around(mnew[index],3)              
    return totalmass

#---Total mass of MW Halo
TM_MWH = ComponentMass(1, "MW_000.txt")
print("the Milky Way mass components")
print("Total Mass of Halo Component:",np.around(TM_MWH,3))
#---Total mass of MW Disk
#TM_MWD = ComponentMass(2, "MW_000.txt")
#print("")
#print("Total Mass of Disk Componet:",np.around(TM_MWD),3)

    
