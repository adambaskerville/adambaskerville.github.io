import numpy as np
from numba import njit
from scipy.special import factorial
import timeit

def cumulative_sum_vectorised_math(a, b, c, d, e, f):
    '''
    This function vectorises the mathematical expression in the nested sums
    '''
    exp_min = 5 - (a + c + d + e)
    exp_max = c + d + f + 1
    exp = 2. ** np.arange(exp_min, exp_max)

    b_fac_grid = np.arange(e)
    g_fac_grid = np.arange(a + b)

    fact_b = factorial(b_fac_grid)
    fact_g = factorial(g_fac_grid)


    return cumulative_sum_JIT(a, b, c, d, e, f, exp_min, exp, fact_b, fact_g)

@njit()
def cumulative_sum_JIT(a, b, c, d, e, f, exp_min, exp, fact_b, fact_g):
    Total = 0
    for ai in range(0, a):
        for bi in range(0, b):
            for ci in range(0, c):
                for di in range(0, d):
                    for ei in range(0, e):
                        for fi in range(0, f):
                            for gi in range(0, ai + bi + 1):
                                for hi in range(0, ci + di + 1):
                                    Total += exp[hi - ei + fi - ai - ci - di + 1 - exp_min] * (ei * ei - 2 * (ei * fi) - 7 * di) * fact_b[bi] * fact_g[gi]
    
    return Total

dim = 5
dim = dim + 1
a = b = c = d = e = f = dim


print(cumulative_sum_vectorised_math(a, b, c, d, e, f))