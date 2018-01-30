"""
Frida Jauregui
ASTR400b HW3 Due Jan 30th
Started Jan 26th

store results in a table and compute total mass of each galaxy
"""

import numpy as np
import astropy.units as u
from astropy.table import Table, Column
from GalaxyMass import ComponentMass

"""
#---notes
store results from GalaxyMass onto a table (no number of particles)
Columms: Galaxy name|Halo Mass[units]|Disk Mass[units]|Bulge[units]|Total Mass|f_bar

Mass in units of 10**12 MSun
Mass of local group(MW, M33 & M31)

Compute: the baryon fraction f_bar = total stellar mass/total(dark+stellar) 
         for each galaxy and the local group
Save: as PDF BONUS if table is created by LaTeX

3.QUESTIONS
"""
tabletxt = 'HW3table.txt'

###Mass Break Down###

#---MW mass components
#---Total mass of MW Halo
TM_MWH = np.around(ComponentMass(1, "MW_000.txt")/1e12, 3)
#---Total mass of MW Disk
TM_MWD = np.around(ComponentMass(2, "MW_000.txt")/1e12, 3)
#---Total mass of MW Bulge
TM_MWB = np.around(ComponentMass(3, "MW_000.txt")/1e12, 3)


#---M31 mass components
#---Total mass of M31 Halo
TM_M31H = np.around(ComponentMass(1, "M31_000.txt")/1e12, 3)
#---Total mass of M31  Disk
TM_M31D = np.around(ComponentMass(2, "M31_000.txt")/1e12, 3)
#---Total mass of M31 Bulge
TM_M31B = np.around(ComponentMass(3, "M31_000.txt")/1e12, 3)


#---M33 mass components
#---Total mass of M33 Halo
TM_M33H = np.around(ComponentMass(1, "M33_000.txt")/1e12, 3)
#---Total mass of M33 Disk
TM_M33D = np.around(ComponentMass(2, "M33_000.txt")/1e12, 3)




#---Total Mass of each Galaxy and Baryon Fraction
TMMW  = np.around(TM_MWH + TM_MWD + TM_MWB, 3)               
TMM31 = np.around(TM_M31H + TM_M31D + TM_M31B, 3)           
TMM33 = np.around(TM_M33H + TM_M33D, 3)
TMLG  = np.around(TMMW + TMM31 + TMM33, 3)
#---Check with print statements
fmw  = np.around((TM_MWD+TM_MWB)/(TMMW), 4)
fm31 = np.around((TM_M31D+TM_M31B)/(TMM31), 4)
fm33 = np.around((TM_M33D)/(TMM33),4)
flg  = np.around((TM_MWD+TM_MWB+TM_M31D+TM_M31B+TM_M33D)/(TMLG), 4)

u = 1e12*u.Msun

#---Table
t = Table()
#---Assign the columns and add them to the table
galacol = Column(name='Galaxy Name', data=['Milky Way','M31','M33'])
halcol  = Column(name='Halo Mass',   data=[TM_MWH,TM_M31H,'0.187'], unit=u)
discol  = Column(name='Disk Mass',   data=[TM_MWD,TM_M31D,'0.009'], unit=u)
bulcol  = Column(name='Bulge Mass',  data=[TM_MWB,TM_M31B,''], unit=u)
totcol  = Column(name='Total Mass',  data=[TMMW,TMM31,'0.196'], unit=u)
t.add_column(galacol)
t.add_column(halcol)
t.add_column(discol)
t.add_column(bulcol)
t.add_column(totcol)

fabrcol = Column(name='f_bar', data=[fmw,fm31,fm33])
t.add_column(fabrcol)

totlgcol =  Column(name='Total Mass of Local Group', data=[TMLG, '', ''], unit=u)
flgcol   =  Column(name='f_bar for Local Group',     data=[flg, '', ''])
t.add_column(totlgcol)
t.add_column(flgcol)

#---Print and save text with LaTex
print(t)
t.write(tabletxt, format='latex', overwrite=True)

"""
#---notes
I get an ValueError: setting an array element with a sequence in M33 vaules
To place the vaules I printed:
M33H, M33D, and TMM33
print("halo mass m33")
print(TM_M33H)
print("disk mass m33")
print(TM_M33D)
print("tot mass m33")
print(TMM33)
and placed them on to the table
"""
