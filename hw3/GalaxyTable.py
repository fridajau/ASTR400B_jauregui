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

###Mass Break Down###

#---MW mass components
#---Total mass of MW Halo
TM_MWH = [np.around(ComponentMass(1, "MW_000.txt"),3)[0]]
##print(TM_MWH)
#---Total mass of MW Disk
TM_MWD = [np.around(ComponentMass(2, "MW_000.txt"),3)[0]]
#---Total mass of MW Bulge
TM_MWB = [np.around(ComponentMass(3, "MW_000.txt"),3)[0]]


#---M31 mass components
#---Total mass of M31 Halo
TM_M31H = np.around(ComponentMass(1, "M31_000.txt"),3)[0]
##print(TM_M31H)
#---Total mass of M31  Disk
TM_M31D = np.around(ComponentMass(2, "M31_000.txt"),3)[0]
#---Total mass of M31 Bulge
TM_M31B = np.around(ComponentMass(3, "M31_000.txt"),3)[0]


#---M33 mass components
#---Total mass of M33 Halo
TM_M33H = [np.around(ComponentMass(1, "M33_000.txt"),3)[0]]
##print(TM_M33H)
#---Total mass of M33 Disk
TM_M33D = [np.around(ComponentMass(2, "M33_000.txt"),3)[0]]

#---Total Mass of each Galaxy and Fbar
TMMW  = TM_MWH + TM_MWD + TM_MWB
TMM31 = TM_M31H + TM_M31D + TM_M31B
TMM33 = TM_M33H + TM_M33D
#TMLG  = TMMW + TMM31 + TM33

#fmw  = (TM_MWD+TM_MWB)/(TMMW) 
fm31 = (TM_M31D+TM_M31B)/(TMM31)
#fm33 = (TM_M33D)/(TM33)
#flg  = (TM_MWD+TM_MWB+TM_M31D+TM_M31B+TM_M33D+TM_M33B)/(TMLG) 


#---Table
t = Table([['Milky Way', 'M31', 'M33'],[TM_MWH,TM_M31H,TM_M33H], [TM_MWD,TM_M31D,TM_M33D]],
          names=('Galaxy Name', 'Halo Mass', 'Disk Mass'))

bulgecol = Column(name='Bulge Mass', data=[TM_MWB,TM_M31B,' '])
t.add_column(bulgecol)

totalmasscol = Column(name='Total Mass', data=[TMMW,TMM31,TMM33])
t.add_column(totalmasscol)

fabrcol = Column(name='f_bar', data=[' ',fm31,' '])
t.add_column(fabrcol)

#t.add_Row(TMLG, 'Total Mass of Local Group')
#t.add_Row(flg, 'f_bar for Local Group')





print(t)


