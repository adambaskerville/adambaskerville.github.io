import numpy as np
import scipy.special

def cumulative_sum_grid_two_sums(a, b, c, d, e, f):
    g = np.arange(a + b).reshape((-1, 1))
    h = np.arange(c + d)

    Total = 0
    for ai in range(0, a):
        for bi in range(0, b):
            gi = g[:ai + bi + 1]
            for ci in range(0, c):
                for di in range(0, d):
                    hi = h[:ci + di + 1]
                    for ei in range(0, e):
                        for fi in range(0, f):
                            Total += np.sum((2.) ** (hi - ei + fi - ai - ci - di + 1) * (ei ** 2 - 2 * (ei * fi) - 7 * di) * scipy.special.factorial(bi) * scipy.special.factorial(gi))
    return Total
dim = 14
dim = dim + 1

cumulative_sum_grid_two_sums(dim, dim, dim, dim, dim, dim)
