import numpy as np
import scipy as sp
import sympy as sy

Z1 = sy.symbols('Z1')
Z2 = sy.symbols('Z2')
Z3 = sy.symbols('Z3')

a1 = sy.symbols('a1')
a2 = sy.symbols('a2')
a3 = sy.symbols('a3')


psi = exp(-a1*z1/2 - a2*z2/2 - a3*z3/2)*F(a1*z1,a2*z2,a3*z3)

T = ep * ( v11*sy.diff(wavefn,z1,z1) + v1*sy.diff(wavefn,z1) + v22*sy.diff(wavefn,z2,z2) + v2*sy.diff(wavefn,z2) + v33*sy.diff(wavefn,z3,z3) + v3*sy.diff(wavefn,z3) + 
v12*sy.diff(wavefn,z1,z2) + v31*sy.diff(wavefn,z3,z1)+v23*sy.diff(wavefn,z2,z3))


V = 2*(2*(Z1*Z2/(z1+z2)+Z1*Z3/(z1+z3)+Z2*Z3/(z2+z3)))*wavefn