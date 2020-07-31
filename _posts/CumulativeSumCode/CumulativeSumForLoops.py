import numpy as np

def cumulative_sum_for_loops(a, b, c, d, e, f):
    '''
    This function calculates the cumulative sum using for loops
    '''
    Total = 0
    for ai in range(0, a):
        for bi in range(0, b):
            for ci in range(0, c):
                for di in range(0, d):
                    for ei in range(0, e):
                        for fi in range(0, f):
                            for gi in range(0, ai+bi+1):
                                for hi in range(0, ci+di+1): 
                                    Total += (2)**(hi-ei+fi-ai-ci-di+1)*(ei**2-2*(ei*fi)-7*di)*np.math.factorial(bi)*np.math.factorial(gi)


    return Total

dim = 14
dim = dim+1

cumulative_sum_for_loops(dim, dim, dim, dim, dim, dim)
