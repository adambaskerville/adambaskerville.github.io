from rogues import *
from matplotlib import pyplot
import seaborn as sns 
from scipy.linalg import eigvals, ordqz
import numpy as np
from numpy import linalg as LA
import mpmath as mp
from mpmath import *

sns.set()
palette = sns.color_palette("bright")

dim = 100
mp.dps = 32

A = lesp(dim)

AT = A.T

Aev  = eigvals(A)
ATev = eigvals(AT)

'''
Aev, Eeg = mp.eig(mp.matrix(A))
ATev, ETeg = mp.eig(mp.matrix(AT))
'''
A_X  = [x.real for x in Aev]
A_Y  = [x.imag for x in Aev]

AT_X = [x.real for x in ATev]
AT_Y = [x.imag for x in ATev]

print(AT_Y)
ax = sns.scatterplot(x=A_X, y=A_Y, color = 'gray', marker='o', label=r'$\mathbf{A}$')
ax = sns.scatterplot(x=AT_X, y=AT_Y, color = 'blue', marker='x', label=r'$\mathbf{A}^T$')

ax.set(xlabel='real', ylabel='imag')
ax.legend()
pyplot.show()