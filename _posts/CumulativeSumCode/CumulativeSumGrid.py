import numpy as np
import scipy.special

def func1(a, b, c, d, e, f):
    # Setup the n-dimensional grid
    ai, bi, ci, di, ei, fi, gi, hi = np.ogrid[:a, :b, :c, :d, :e, :f, :a + b - 1, :c + d - 1]
    # Calculate the mathematics within the summations
    Total = (2.) ** (hi - ei + fi - ai - ci - di + 1) * (ei ** 2 - 2 * (ei * fi) - 7 * di) * scipy.special.factorial(bi) * scipy.special.factorial(gi)
    # Mask out of range elements for last two inner loops
    mask = (gi < ai + bi + 1) & (hi < ci + di + 1)
    return np.sum(Total * mask)

dim = 9
dim = dim + 1
func1(dim, dim, dim, dim, dim, dim)
