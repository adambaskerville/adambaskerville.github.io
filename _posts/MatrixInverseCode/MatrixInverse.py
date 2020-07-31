import scipy.linalg as spla
import numpy as np
M11 = 1.01
M12 = 1.00 
M13 = 1.00

M21 = 1.00
M22 = 1.01
M23 = 1.00

M31 = 1.00
M32 = 1.00
M33 = 1.00

A = np.array([[M11, M12, M13], 
              [M21, M22, M23],
              [M31, M32, M33]])

b = np.array([[4], [7.9999999999999999]])

def np_inv(A, b):
    return np.linalg.inv(A)

def np_solve(A, b):
    return np.linalg.solve(A, b)

def svd_inv(A, b):
    #u, s, v = np.linalg.svd(A)
    #Ainv = np.dot(v.transpose(), np.dot(np.diag(s**-1), u.transpose()))
    Ainv = np.linalg.pinv(A)
    return Ainv

def svd_solve(A, b):
    U, S, VT = np.linalg.svd(A)

    C = np.dot(U.T, b)

    w = np.linalg.solve(np.diag(S), C)

    x = VT.T @ w

    return x

print("np_inv\n", np_inv(A, b))

print("svd_inv\n", svd_inv(A, b))

#print("np_solve\n", np_solve(A, b))

#print("SVD_solve\n", svd_solve(A, b))