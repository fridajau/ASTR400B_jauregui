import numpy as np
import astropy.units as u
import matplotlib
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterofMass import CenterOfMass

#---class
#Create class that will take a galaxy, snapshot, and particle type
class SolarParticles:
    def __init__(self, galaxy, snap, ptype):
        #--add string to the text file format "000"
        #--remove all but the last 3 digits
        ilbl = '000' + str(snap)
        ilbl = ilbl[-3:]
        #--create filenames
        self.filename='%s_'%(galaxy) + ilbl + '.txt'
        #--use HighRes
        #--creating a COM for M31 using Disk Particles from CenterOfMass
        self.COM_M31 = CenterOfMass('HighRes/'+self.filename,2)
        M31_pos = self.COM_M31.COM_P(1.0,4.0)
        M31_vel = self.COM_M31.COM_V(M31_pos[0],M31_pos[1],M31_pos[2])
        
        self.COMX = self.COM_M31.x - float(M31_pos[0]/u.kpc)
        self.COMY = self.COM_M31.y - float(M31_pos[1]/u.kpc)
        self.COMZ = self.COM_M31.z - float(M31_pos[2]/u.kpc)
        
        self.COMVX = self.COM_M31.vx - float(M31_vel[0]/(u.km/u.s))
        self.COMVY = self.COM_M31.vy - float(M31_vel[1]/(u.km/u.s))
        self.COMVZ = self.COM_M31.vz - float(M31_vel[2]/(u.km/u.s))
        #--COM object at snapshot 0
        COM_M310 = CenterOfMass('HighRes/M31_000.txt',2)
        M31_pos0 = COM_M310.COM_P(1.0,4.0)
        M31_vel0 = COM_M310.COM_V(M31_pos0[0],M31_pos0[1],M31_pos0[2])
        #--x,y,z pos & vel of m31 rel to COM
        COMX0 = COM_M310.x - float(M31_pos0[0]/u.kpc)
        COMY0 = COM_M310.y - float(M31_pos0[1]/u.kpc)
        COMZ0 = COM_M310.z - float(M31_pos0[2]/u.kpc)
        #--radial distance from center of galaxy
        RadPos = np.sqrt(COMX0**2 + COMY0**2)
        #--index the vaules within 7-9 kpc CHANGE ZCOM to a resaonable scale height
        self.Rindex = np.where((RadPos > 7) & (RadPos < 9) & (COMZ0 > -5) & (COMZ0 < 5))



        
        #--functions
    #Obtain radial positions not at com position for any galaxy
    def RadialPos(self):
        #INPUT:   self
        #RETURNS: pos with the radial index
        
        return self.COM_M31.x[self.Rindex], self.COM_M31.y[self.Rindex], self.COM_M31.z[self.Rindex]

    
    #Obtain radial positions AT com pos & vel(if needed)
    def RadialPosCOM(self):
        #INPUT:   self  
        #RETURNS: 3d coord of pcom or vcom of particals within 7-9 kpc

        #--use the Rindex to always return the ring of particles at COM of x,y,z,vx,vy,vz
        nM31x = self.COMX[self.Rindex]
        nM31y = self.COMY[self.Rindex]
        nM31z = self.COMZ[self.Rindex]

        nM31vx = self.COMVX[self.Rindex]
        nM31vy = self.COMVY[self.Rindex]
        nM31vz = self.COMVZ[self.Rindex]

        return nM31x, nM31y, nM31z, nM31vx, nM31vy, nM31vz

    
    #Obtain a vaule of radii from M33 that will pass through/bound the ring particals as MW-M31 merges
    def IntoM33(self):
        #INPUT:   self
        #RETURNS: seperation of x - m33com.x of pos vaules

        #--use center of mass to find half the radius of m33
        #--find m33 com from center of mass
        #--use radpos to get the x,y,z positions of m33
        #--rindex for range of radii
        #--delta_r must be less than r/2 of m33
        #--returns some seperation of x - m33com.x
    
    

 """   
#Obtain positions of particles for M31, MW, M33
#--m31
Ring_M31    = SolarParticles("M31", 700, 2)  #--"galaxy", snapshot, disk particles
radial_ring = Ring_M31.RadialPosCOM()        #--RadialCOM() centering M31
xm31  = np.around(radial_ring[0],3)
ym31  = np.around(radial_ring[1],3)
zm31  = np.around(radial_ring[2],3)

#--m31 all
Disk_M31= CenterOfMass("M31_000.txt", 2)
M31_pos = Disk_M31.COM_P(1.0,4.0)
Dxm31 = Disk_M31.x #- float(MW_pos[0]/u.kpc)   #--commented out COM for real locations at given snapshot 
Dym31 = Disk_M31.y #- float(MW_pos[1]/u.kpc)    
Dzm31 = Disk_M31.z #- float(MW_pos[2]/u.kpc)

#--mw
Disk_MW = CenterOfMass("MW_000.txt", 2)
MW_pos  = Disk_MW.COM_P(1.0,4.0)
xmw = Disk_MW.x #- float(MW_pos[0]/u.kpc)
ymw = Disk_MW.y #- float(MW_pos[1]/u.kpc)
zmw = Disk_MW.z #- float(MW_pos[2]/u.kpc)

#--m33
Disk_M33 = CenterOfMass("M33_000.txt", 2)
M33_pos  = Disk_M33.COM_P(1.0,4.0)
xm33 = Disk_M33.x #- float(M33_pos[0]/u.kpc)
ym33 = Disk_M33.y #- float(M33_pos[1]/u.kpc)
zm33 = Disk_M33.z #- float(M33_pos[2]/u.kpc)


#--plots
#Make a histogram for M31 vs MW (look at M33)
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

#--hist for M31
plt.hist2d(xm31, ym31, bins=600, norm=LogNorm(), cmap='magma')
plt.colorbar()
#--hist for MW
#plt.hist2d(xmw, ymw, bins=650, norm=LogNorm(), cmap='Reds')
#plt.colorbar()
#--hist for M33
#plt.hist2d(xm33, ym33, bins=650, norm=LogNorm(), cmap='Blues')
#plt.colorbar()
#--label for snap shot
#plt.text(-50, 60,'T=3.87')
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('y (kpc)', fontsize=22)
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size
#plt.ylim(-150,150)
#plt.xlim(-130,130)
plt.show()

#########################################################################################

####Ratio of Sun Candidates at a Given Snapshot####

#--selected particles at a choosen snapshot, from above initialization(700), at COM
ringCOM   =  Ring_M31.RadialPosCOM()
ring_part = np.sqrt(ringCOM[0]**2 + ringCOM[1]**2 + ringCOM[2]**2)
#--all ring particles
ring_tot  = np.size(ringCOM[0])
print("magnitude of ring of part",ring_tot)

#--initialize a new object at snap0 COM
diskm31_S = SolarParticles("M31", 000, 2)
ringcom_S    = diskm31_S.RadialPosCOM()
ringofpart_S = np.sqrt(ringcom_S[0]**2 + ringcom_S[1]**2 + ringcom_S[2]**2)

#--plots
#Make a histogram plot thatll show the radial distribution of candidate suns w/
#respect to the center of the MW & M31 remnant
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
plt.xlabel('Distance rom Merger Remnant (kpc)', fontsize=22)
plt.ylabel('Solar Candinates', fontsize=22)

#--hist of ring of particles
plt.hist(ring_part, bins=120, color='black', histtype="step")
#--hist of disk particles at snapshot0 
plt.hist(ringofpart_S, bins=120, color='blue', histtype="step")
#--add lines to label plot
plt.axvline(8.29, color='red', ls='--')
plt.text(8.29,250,'Current Distance [8.29 kpc]',rotation=90)
plt.text(40,90,'Black: T=10.0')
plt.text(40,80,'Blue:  T=0.0')
#plt.xticks([20,40,60,80,100])
#plt.axis([0, 150, 0, 400])
plt.show()
"""
