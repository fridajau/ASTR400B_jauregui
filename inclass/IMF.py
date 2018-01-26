"""
Frida Jauregui
ASTR400b
In class worksheet 1
Jan 25, 2018
"""
from scipy.integrate import quad
import numpy as np

"""
#---goal
complete the IMF of Salpeter
"""

#---IMF
def Salpeter(M, Mmin, Mmax):
    Z = quad(lambda M: M**(- 2.35), Mmin, Mmax)  #determine the magnitude of the integral
    A = 1/Z[0]                                   #normalize to 1
    return A*M**-2.35                            #return the normalized function


#---determine the fraction of stars that are greater than one to 120 solar masses
Frac = quad(lambda M: Salpeter(M, 0.1, 120), 1.0, 120)
#---want the first number as the second one is an error vaule
print("Fractional number of Stars more masssive than our Sun:", np.around(Frac[0], 3))
