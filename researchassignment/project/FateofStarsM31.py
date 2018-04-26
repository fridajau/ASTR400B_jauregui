import numpy as np
import astropy.units as u
import matplotlib
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from ReadFile import Read
from CenterofMass import CenterOfMass
from Orbits import OrbitCOM

#---class
#---index disk particles from data to work only with disk particles
class SolarParticles:
    def __init__(self, filename, ptype):
        self.time, self.total, self.data = Read(filename)
        #--storing data for disk particles
        self.index = np.where(self.data['type'] == ptype)
        
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]

#--Index data to return vaules of M31 within a radial distance of 7-9 kpc.
#--Contain  particles that are within half the scale hieght of M31
#--Get COM of the orbits of the galaxies and change the particle pos/vel vector to the COM frame for
#each snap shot <--obtain files from this
#--Compute the circular velocity of M31 at each pericenter and determine what percentage of stars will get
#transfered or pass through the milkyway <-- will need seperate data sets to tell the difference b/t particles from
#                                            the galaxies, maybe visual aid
#                                            if time allows M33 as well

#--Use lab7 to help rearrange the merger to be face-on and plot the density profile at snap shots
#--Create a histogram to see the radial distribution with repect to M31's galatic center of solar particles

    def RadialIndex(self, x,y,z):
        #INPUT:  x,y,z postions
        #OUTPUT: the ring of particles

        #--use CenterOfMass to obtain the COM of x,y,z vx,xy,vz of M31
        #--snap shots to look out for (0.0, 3.87, 5.87, 6.2, 10.0)
        #--obtain the index ring of particles to use later on to define at wanted snap shots
        #--follow this function for radial positions

        #--creating a COM obeject for M31 using Disk Particles from CenterOfMass
        COM_M31 = CenterOfMass("M31_000.txt",2)
        M31_pos = COM_M31.COM_P(1.0,4.0)
        M31_vel = COM_M31.COM_V(M31_pos[0],M31_pos[1],M31_pos[2])
        #--x,y,z pos & vel of m31 rel to COM
        COMX = COM_M31.x - float(M31_pos[0]/u.kpc)
        COMY = COM_M31.y - float(M31_pos[1]/u.kpc)
        COMZ = COM_M31.z - float(M31_pos[2]/u.kpc)
        
        COMVX = COM_M31.vx - float(M31_vel[0]/(u.km/u.s))
        COMVY = COM_M31.vy - float(M31_vel[1]/(u.km/u.s))
        COMVZ = COM_M31.vz - float(M31_vel[2]/(u.km/u.s))
        #--3d radial position from galatic center
        RadPos = np.sqrt(COMX**2 + COMY**2)
        #--index the vaules within 7-9 kpc
        Rindex = np.where((RadPos > 7) & (RadPos < 9) & (COMZ > -1.0) & (COMZ < 1.0))
    
        return Rindex
    
    def RadialPos(self, x,y,z):
        #INPUT: x,y,z positions, center of galaxy=0  
        #RETURNS: 3d coord of pcom or vcom of particals within 7-9 kpc a snap shot 0

        #--index M31 vaules and return an array with the desired particles

        COM_M31 = CenterOfMass("M31_000.txt",2)
        M31_pos = COM_M31.COM_P(1.0,4.0)
        M31_vel = COM_M31.COM_V(M31_pos[0],M31_pos[1],M31_pos[2])
       
        COMX = COM_M31.x - float(M31_pos[0]/u.kpc)
        COMY = COM_M31.y - float(M31_pos[1]/u.kpc)
        COMZ = COM_M31.z - float(M31_pos[2]/u.kpc)
        
        COMVX = COM_M31.vx - float(M31_vel[0]/(u.km/u.s))
        COMVY = COM_M31.vy - float(M31_vel[1]/(u.km/u.s))
        COMVZ = COM_M31.vz - float(M31_vel[2]/(u.km/u.s))

        #--3d radial position from galatic center
        RadPos = np.sqrt(COMX**2 + COMY**2)
        #--use Rindex function to always return the ring of particles
        Rindex = self.RadialIndex(COMX,COMY,COMZ)
        nM31x = COMX[Rindex]
        nM31y = COMY[Rindex]
        nM31z = COMZ[Rindex]

        return nM31x, nM31y, nM31z
    
#--checking if it returns ring of particles
Disk_M31= SolarParticles("M31_000.txt", 2)
radial_ring = Disk_M31.RadialPos(-377,608,-284)
testx = radial_ring[0]
testy = radial_ring[1]
testz = radial_ring[2]

rad = np.average(np.sqrt(testx**2 + testy**2 + testz**2))
print("ave pos", rad)

fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)

plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('y (kpc)', fontsize=22)

plt.hist2d(testx, testy, bins=400, norm=LogNorm(), cmap='magma')
plt.colorbar()
plt.ylim(-30,30)
plt.xlim(-30,30)
plt.show()

"""
    def M31SunsOrbit(self, DiskPart, Snap, end, n):
        #INPUT:  the particles I'm going to orbit, snap shots, integrator
        #OUTPUT: return the orbit of the M31 & MW system at deseried snap shots
            
        #--Compute the orbit of my ring particles with MW
        #--I have M31_Orbit.txt which contains the orbit of disk particles at all snap shots
        #--Use Orbits.py as a template
        #--Snaps shots to obtain (0.0, 3.87, 5.87, 6.2 & 10.0)
        #--file of the orbit of MW and the reduced M31 particles
        #--take steps from Orbits.py and hw5
        
        #--define the filename for the file that will store the orbit
        fileout = "Orbit_{}.txt".format(DiskPart)

        #--define an array with 7 columns that will store
        #--t, x,y,z of the COM
        Orbit = np.zeros((int(end/n)+1, 4))

        #--for loop from start to end+1 in intervals of n
        for i in np.arange(Snap, end+n, n):
            #--define the filename for the galaxy using hw5 as example
            ilbl = '000' + str(i)
            ilbl = ilbl[-3:]
            filename = "%s_"%(DiskPart) + ilbl + '.txt'
            
            #--COM M31 Sun-like stars
            M31 = np.genfromtxt('M31_lifetime.txt', dtype = float, names=True)
            COM = CenterOfMass('VLowRes/'+filename, 2)
            Disk_M31 = SolarParticles("M31_000.txt", 2)
            D_COMM31 = Disk_M31.RadialPos(-377, 608, -284)

            #--row index of the Orbit array given as int(i)/n)
            Orbit[int(i/n),0] = float(COM.time/u.Myr)/1000
            #--store the COM pos and vel in the Orbit array, divide out the units
            Orbit[int(i/n),1] = float(D_COMM31[0])
            Orbit[int(i/n),2] = float(D_COMM31[1])
            Orbit[int(i/n),3] = float(D_COMM31[2])
            
        fileout = 'M31_Ring_0.txt'
    
        #--save the array Orbit to a file
        np.savetxt(fileout, Orbit, header='t, x, y, z', comments='#',
                   fmt=['%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f', '%.2f'])

        return Orbit
    

    #def SunCand_M33(self, ):
        #--find the percentage of Sun Candinate stars that will fall into M33 from the Merger (3.87, 5.87, 6.2, 10.0)
        #--use/find the radial vel at 8 kpc of M31
        #--
"""     


    

        
