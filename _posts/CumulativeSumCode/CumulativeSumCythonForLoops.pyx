#cython: language_level=3
cimport cython
import numpy as np

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cpdef cumulative_sum_for_loops(int a, int b, int c, int d, int e, int f):
    '''
    This function calculates the cumulative sum using for loops
    '''
    cdef int ai, bi, ci, di, fi, gi, hi
    cdef long double Total

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