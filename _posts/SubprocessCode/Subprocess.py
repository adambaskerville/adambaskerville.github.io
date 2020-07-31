import os
import signal
from subprocess import Popen, PIPE, STDOUT
import scipy as sp
import numpy as np
from scipy.optimize import minimize
      
def CallCpp(params):
    x = params[0]
    y = params[1]

    # Pipe the x and y values to the C++ executable
    process.stdin.write('{} {}\n'.format(x, y))
    process.stdin.flush()
    
    FuncValue = float(process.stdout.readline())

    print("x={:.6E} y={:.6E} f(x,y)={:.16f}\n".format(x, y, FuncValue))

    return FuncValue

# Initiate the subprocess by calling the Subprocess executable we just compiled
process = Popen('./Subprocess', stdin=PIPE, stdout=PIPE, universal_newlines=True, shell=True, preexec_fn=os.setsid)

# Tell Python that the first line output by C++ is: "This C++ program will be called from Python continuously" otherwise it will prevent the values of x and y being piped
print(process.stdout.readline()) 

# provide the starting values for x and y
initial_guess = [2,2]

# Minimise the functin using the Powell algorithm
Func_min = minimize(CallCpp, initial_guess, method='nelder-mead', options={'ftol': 1e-10,'disp': True})

# Kill the subprocess
os.killpg(os.getpgid(process.pid), signal.SIGTERM)