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

M33 doesn't have a bulge component
"""
def ComponentMass(PType, filename):
    #---to be read in the file
    time, total, data = Read(filename)

    #---an array to store Ptype
    index = np.where(data['type'] == PType)

    #---Mass
    mnew = data['m'][index]*1e10*u.Msun
    totalmass = np.around(mnew,3)              
    return totalmass
"""
#---check work
print("The Milky Way mass components")
print("")

#---Total mass of MW Halo
TM_MWH = ComponentMass(1, "MW_000.txt")
print("Total Mass of Halo Component:",np.around(TM_MWH[0],3))

#---Total mass of MW Disk
TM_MWD = ComponentMass(2, "MW_000.txt")
print("Total Mass of Disk Componet:",np.around(TM_MWD[0],3))

#---Total mass of MW Bulge
TM_MWB = ComponentMass(3, "MW_000.txt")
print("Total Mass of Bulge Component:", np.around(TM_MWB[0],3))

print("")
print("M31 mass components")
print("")

#---Total mass of M31 Halo
TM_M31H = ComponentMass(1, "M31_000.txt")
print("Total Mass of Halo Component:",np.around(TM_M31H[0],3))

#---Total mass of M31  Disk
TM_M31D = ComponentMass(2, "M31_000.txt")
print("Total Mass of Disk Component:",np.around(TM_M31D[0],3))

#---Total mass of M31 Bulge
TM_M31B = ComponentMass(3, "M31_000.txt")
print("Total Mass of Bulge Component:",np.around(TM_M31B[0],3))

print("")
print("M33 mass components")

#---Total mass of M33 Halo
TM_M33H = ComponentMass(1, "M33_000.txt")
print("Total Mass of Halo Component:",np.around(TM_M33H[0],3))

#---Total mass of M33 Disk
TM_M33D = ComponentMass(2, "M33_000.txt")
print("Total Mass of Disk Component:",np.around(TM_M33D[0],3))
"""
