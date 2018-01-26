"""
Frida Jauregui
ASTR400b
In class worksheet 1
Jan 25, 2018
"""
import numpy as np
import astropy.units as u
from pylab import *


"""
#---goal
plot the schechter function
"""

def Schechter(M, pstar, Mstar, alpha):
    phi = (0.4*np.log(10))*pstar*10**(0.4*(Mstar - M)*(alpha + 1))*exp(-10**(0.4*(Mstar - M)))
    return phi

#---constants
pstar = 1.66*10**(-2)
Mstar = -23.19

#---plot
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)
MK = np.arange(-26, -17, 0.1)                       #an array to store K-band Mag                    
plt.semilogy(MK,Schechter(MK,pstar,Mstar,-0.81), color='blue', linewidth=5, label='Smith+09')

#---plot alphas
plt.semilogy(MK,Schechter(MK,pstar,Mstar,-0.6), color='blue', linestyle=":", linewidth=3, label=r'low $\alpha$')
plt.semilogy(MK,Schechter(MK,pstar,Mstar,-1.35), color='blue', linestyle="--",linewidth=3, label=r'high $\alpha$')

#---plot labels
plt.xlabel(r'M$_k$ + 5Log($h$)', fontsize=22)
plt.ylabel(r'$\Phi$ (Mpc$^{-3}h^3$/mag)', fontsize=22)

#---limits
plt.xlim(-17,-26)

#---lengend
legend = ax.legend(loc='upper right',fontsize='x-large')

plt.show()

# Save to a file
ax.set_rasterized(True)
plt.savefig('Schechter_Jan25.eps', rasterized=True, dpi=350)




