from rogues import *
from matplotlib import pyplot
import seaborn as sns 
from scipy.linalg import eigvals, ordqz
import numpy as np
from numpy import linalg as LA
from flint import *

sns.set()
palette = sns.color_palette("bright")
dim = 100
# get Grcar matrix
A = grcar(dim)
AT = A.T
B = np.identity(dim)
BT = B.T
#Aev = eigvals(A)
#ATev = eigvals(A.T)

AA, BB, alpha, beta, Q, Z = ordqz(A, B)

Aev = np.divide(alpha, beta)

AAT, BBT, alphaT, betaT, QT, ZT = ordqz(AT, BT)

ATev = np.divide(alphaT, betaT)

A_X = [x.real for x in Aev]
A_Y = [x.imag for x in Aev]

AT_X = [x.real for x in ATev]
AT_Y = [x.imag for x in ATev]


# compute pseudospectrum for the levels of interest between [1e-5, 1]
#pseudo = NonnormalAuto(A, 1e-5, 1)

# plot
#pseudo.plot([10**k for k in range(-4, 0)], spectrum=ANS2)
ax = sns.scatterplot(x=A_X, y=A_Y, color = 'gray', marker='o', label=r'$\mathbf{A}$')
ax = sns.scatterplot(x=AT_X, y=AT_Y, color = 'blue', marker='x', label=r'$\mathbf{A}^T$')

ax.set(xlabel='real', ylabel='imag')
ax.legend()

pyplot.show()
'''
pyplot.scatter(X,Y, color='red')
pyplot.scatter(X2,Y2, color='blue')
pyplot.show()
'''