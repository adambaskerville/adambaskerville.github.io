from scipy.optimize import minimize
import numpy as np
def min(c):
    E = (3/2)*c - (2*(2**0.5)*(c**0.5))/(np.pi**0.5)
    
    #E = (c**2/2) - c
    print("c = {:.6f}, Energy = {:.16f}".format(c[0], E[0]))
    return E

# Provide an initial guess for the variational parameter c
c = 1

bnds = ((0,10),)
# Optimise the min function by varying c with a tolerance of 1e-6
hydrogen_gs = minimize(min, c, tol=1e-6, bounds=bnds) 