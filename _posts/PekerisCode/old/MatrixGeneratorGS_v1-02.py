import os
import cexprtk
import numpy as np
from numpy import genfromtxtimport multiprocessing
import parmap
import scipy
from scipy.linalg import eigvalsh, ordqz
import time

from sympy import Symbol, solve
np.seterr(divide='ignore')
# Specify number of workers for multiprocessing
number_of_workers = multiprocessing.cpu_count()+1

Global_dict = {"A"  : 1.1,
               "B"  : 1.1,
               "C"  : 2.2,
               "K"  : 1,
               "ep" : 1,
               "t"  : 1,
               "p"  : 1,
               "l"  : 1,
               "m"  : -1,
               "n"  : 0,
               "lp" : 1,
               "mp" : 1,
               "np" : 1,
               "m1" : 1,
               "m2" : 1,
               "m3" : 10000,
               "E"  : 1,
               "Z1" : -1,
               "Z2" : -1,
               "Z3" : 2,
               "ss" : 1,
               "hh" : 1,
               "1" : 1,
               "AllEigVals" : False,
               "Root" : 2,
               "EigenProb" : "Energy"
               }

# Setup the cexprtk symbol table
st = cexprtk.Symbol_Table({"A"  : Global_dict["A"],
                           "B"  : Global_dict["B"],
                           "C"  : Global_dict["C"],
                           "K"  : Global_dict["K"],
                           "ep" : 1,
                           "t"  : 1,
                           "p"  : 1,
                           "l"  : Global_dict["l"],
                           "m"  : Global_dict["m"],
                           "n"  : Global_dict["n"],
                           "lp" : Global_dict["lp"],
                           "mp" : Global_dict["mp"],
                           "np" : Global_dict["np"],
                           "m1" : Global_dict["m1"],
                           "m2" : Global_dict["m2"],
                           "m3" : Global_dict["m3"],
                           "E"  : 1,
                           "Z1" : Global_dict["Z1"],
                           "Z2" : Global_dict["Z2"],
                           "Z3" : Global_dict["Z3"],
                           "ss" : 1,
                           "hh" : 1,
                           }, 
                           add_constants=True)

# Read in the recursion relations
RRHH = genfromtxt(os.path.expanduser('~/Postdoc/post_doc_2019/3Body2/3Body2FC/RR_HH.txt'), dtype=(int,int,int,'S1000'))
RRSS = genfromtxt(os.path.expanduser('~/Postdoc/post_doc_2019/3Body2/3Body2FC/RR_SS.txt'), dtype=(int,int,int,'S1000'))

# Store the recursion relations in a dictionary and compile them into run time code using cexprtk
RRDictHH = {(row[0],row[1],row[2]) : cexprtk.Expression(str(row[3], 'utf-8'), st) for row in RRHH}
RRDictSS = {(row[0],row[1],row[2]) : cexprtk.Expression(str(row[3], 'utf-8'), st) for row in RRSS}

mat_size = 2856

def testlmn(l, m, n, A, B, C, m1, m2, m3, Z1, Z2, Z3):
    return 1/4*(-6*((1/3+n+n**2)*m3*(m1+m2)*B**2+4/3*C*(n+1/2)*m3*(m+1/2)*m1*B+C**2*(m2+m3)*(m+m**2+1/3)*m1)*(l+1/2)*A**3+(-6*(1/3+n+n**2)*m3*(m1+m2)*(m+1/2)*B**3-12*m3*(1/6*(n+1/2)*(m1+m2)*(l**2+m**2+l+m)*C+Z1*Z2*m1*m2*(1/3+n+n**2))*B**2-8*(1/4*(l**2+n**2+l+n)*(m2+m3)*C+(n+1/2)*((Z2+Z3)*Z1+Z2*Z3)*m3*m2)*C*m1*(m+1/2)*B-12*(m+m**2+1/3)*(1/2*(n+1/2)*(m2+m3)*C+m2*m3*Z1*Z3)*C**2*m1)*A**2-8*((n+1/2)*m3*(m+1/2)*B**2+(1/4*(m**2+n**2+m+n)*(m1+m3)*C+(n+1/2)*m1*((Z2+Z3)*Z1+Z2*Z3)*m3)*B+((n+1/2)*C+((Z2+Z3)*Z1+Z2*Z3)*m3)*C*m1*(m+1/2))*C*(l+1/2)*m2*B*A-12*C**2*m2*(1/2*(m+1/2)*(m1+m3)*B+1/2*(n+1/2)*(m1+m3)*C+m1*m3*Z2*Z3)*B**2*(l**2+l+1/3))/B**3/C**3/A**3/m1/m2/m3


l =  1
m =  -1
n =  0
A =  1.1
B =  1.1
C =  2.2
m1 = 1.0
m2 = 1.0
m3 = 10000.0
Z1 = -1.0
Z2 = -1.0
Z3 = 2.0

start = time.time()
for i in range(100000):
    RRCython.RRHH(l, m, n, m1, m2, m3, Z1, Z2, Z3, A, B, C)

end = time.time()
print("Took {} seconds.".format(end - start))

start = time.time()
for i in range(100000):
    RRDictHH[1,-1,0]()

end = time.time()
print("Took {} seconds.".format(end - start))

#print(testlmn(l, m, n, A, B, C, m1, m2, m3, Z1, Z2, Z3))

exit()
## Symmetric numbering scheme
size = 1.4*mat_size
omega = Symbol('omega', real=True)
omega = int(round(solve(15/16+(17/12)*omega+(5/8)*omega**2+(1/12)*omega**3+(1/16)*(1)**omega-size, omega)[0]))

LMN_MAP_SYM = []    
for ww in range(omega + 1):
    for vv in range(ww + 1):
        for uu in range(vv + 1):
            l = uu
            m = vv - uu
            n = ww - vv
            if l+m+n <= omega and l<=m:
                lmn = (l,m,n) 
                LMN_MAP_SYM.append(lmn) 

omega = ((1/3)*(81*size+3*(729*size**2-3)**0.5)**(1/3)+1/(81*size+3*(729*size**2-3)**0.5)**(1/3)-2)
omega = int(round(omega, 0))

LMN_MAP_ASYM = []    
for ww in range(omega + 1):
    for vv in range(ww + 1):
        for uu in range(vv + 1):
            l = uu
            m = vv - uu
            n = ww - vv
            if l + m + n <= omega:
                lmn = (l,m,n) 
                LMN_MAP_ASYM.append(lmn) 

## Antisymmetric numbering scheme

omega = Symbol('omega', real=True)
size = 2*mat_size
omega = int(round(solve(15/16+(17/12)*(omega-1)+(5/8)*(omega-1)**2+(1/12)*(omega-1)**3+(1/16)*(1)**(omega-1)-size, omega)[0]))
LMN_MAP_ANTISYM = []    
for ww in range(omega+1):
    for vv in range(ww+1):
        for uu in range(vv+1):
            l = uu
            m = vv - uu
            n = ww - vv
            if l+m+n <= omega and l<m:
                lmn = (l,m,n) 
                LMN_MAP_ANTISYM.append(lmn) 

def mat_build_ASYM(H, mat_size, rrtype):
    '''

    '''
    if rrtype == 'HH':
        dictionary = RRDictHH
    elif rrtype == 'SS':
        dictionary = RRDictSS

    mat_elem_list     = [0] * (mat_size)

    st.variables['l'] = LMN_MAP_ASYM[H][0]
    st.variables['m'] = LMN_MAP_ASYM[H][1]
    st.variables['n'] = LMN_MAP_ASYM[H][2]

    for J in range(H+1):
        lam = LMN_MAP_ASYM[J][0] - LMN_MAP_ASYM[H][0]
        mu  = LMN_MAP_ASYM[J][1] - LMN_MAP_ASYM[H][1]
        nu  = LMN_MAP_ASYM[J][2] - LMN_MAP_ASYM[H][2]

        if (lam, mu, nu) not in dictionary:
            mat_elem_list[J] += 0
        else:
            mat_elem_list[J] += dictionary[lam, mu, nu]()

    return np.array(mat_elem_list, dtype='longdouble').reshape(mat_size,1)

def mat_build_SYM(H, mat_size, rrtype):
    if rrtype == 'HH':
        dictionary = RRDictHH
    elif rrtype == 'SS':
        dictionary = RRDictSS

    mat_elem_list = [0] * (mat_size)
    
    for J in range(H, mat_size):
        # Specify what L, M and N are to the cexprtk symbol table
        st.variables['l'] = L =  LMN_MAP_SYM[H][0]
        st.variables['m'] = M =  LMN_MAP_SYM[H][1]
        st.variables['n'] = N =  LMN_MAP_SYM[H][2]
        # Specify what Lp, Mp and Np are
        Lp =  LMN_MAP_SYM[J][0]
        Mp =  LMN_MAP_SYM[J][1]
        Np =  LMN_MAP_SYM[J][2]
        # Take their difference
        lam = Lp - L 
        mu  = Mp - M
        nu  = Np - N
        # Apply symmetry rules
        if (L == M) and (Lp == Mp):
            if (lam, mu, nu) not in dictionary:
                newval1 = 0
            else:
                newval1 = dictionary[lam, mu, nu]()/2
        else:
            if (lam, mu, nu) not in dictionary:
                newval1 = 0
            else:
                newval1 = dictionary[lam, mu, nu]()
        # Apply symmetry rules
        if (L != M) and (Lp != Mp):
            lam = Lp - M
            mu  = Mp - L
            nu  = Np - N
            # Update values in symbol table
            st.variables['l'] = M 
            st.variables['m'] = L 
            st.variables['n'] = N 

            if (lam, mu, nu) not in dictionary:
                newval2 = 0
            else:
                newval2 = dictionary[lam, mu, nu]()
        else:
            newval2 = 0

        # Add to list. it is a factor of 8 out due to the volume element. This is included here, but is not needed as just scales the matrix which does not effect the eigenvalues
        mat_elem_list[J] += newval1 + newval2
        
    return np.array(mat_elem_list, dtype='longdouble').reshape(mat_size,1)

def mat_build_ANTISYM(H, mat_size, rrtype):
    if rrtype == 'HH':
        dictionary = RRDictHH
    elif rrtype == 'SS':
        dictionary = RRDictSS

    mat_elem_list = [0] * (mat_size)

    for J in range(H, mat_size):
        # Specify what L, M and N are to the cexprtk symbol table
        st.variables['l'] = L =  LMN_MAP_ANTISYM[H][0]
        st.variables['m'] = M =  LMN_MAP_ANTISYM[H][1]
        st.variables['n'] = N =  LMN_MAP_ANTISYM[H][2]
        # Specify what Lp, Mp and Np are
        Lp =  LMN_MAP_ANTISYM[J][0]
        Mp =  LMN_MAP_ANTISYM[J][1]
        Np =  LMN_MAP_ANTISYM[J][2]
        # Tabke their difference
        lam = Lp - L 
        mu  = Mp - M
        nu  = Np - N
        # Apply symmetry rules
        if (L == M) and (Lp == Mp):
            if (lam, mu, nu) not in dictionary:
                newval1 = 0
            else:
                newval1 = dictionary[lam, mu, nu]()/2
        else:
            if (lam, mu, nu) not in dictionary:
                newval1 = 0
            else:
                newval1 = dictionary[lam, mu, nu]()
        # Apply symmetry rules
        if (L != M) and (Lp != Mp):
            lam = Lp - M
            mu  = Mp - L
            nu  = Np - N
            # Update values in symbol table
            st.variables['l'] = M 
            st.variables['m'] = L 
            st.variables['n'] = N 

            if (lam, mu, nu) not in dictionary:
                newval2 = 0
            else:
                newval2 = dictionary[lam, mu, nu]()
        else:
            newval2 = 0
            
        # Add to list. it is a factor of 8 out due to the volume element. This is included here, but is not needed as just scales the matrix which does not effect the eigenvalues
        mat_elem_list[J] += 8*(newval1 - newval2)
        
    return np.array(mat_elem_list, dtype='longdouble').reshape(mat_size,1)

# Build a list to be used for parallelization in parmap
Jrange = range(mat_size) # For parmap parallelization
Jlist = [*Jrange]

start = time.time()

listHH = parmap.map(mat_build_SYM, Jlist, mat_size, 'HH', pm_processes=number_of_workers)
listSS = parmap.map(mat_build_SYM, Jlist, mat_size, 'SS', pm_processes=number_of_workers)

HHMat = np.hstack(listHH)
SSMat = np.hstack(listSS)
#print(HHMat)

end = time.time()
print("Building of the {} x {} Hamiltonian and Overlap matrices took {} seconds.".format(mat_size, mat_size, end - start))
ev = -eigvalsh(HHMat, (-1)*SSMat)[-1]
print(ev)
