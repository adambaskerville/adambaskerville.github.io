---
layout: post
title: "T>T: Call C++ From Python Without Wrapping"
date: 2019-08-29
excerpt: "Using subprocess to quickly and easily call a C++ program, whilst keeping input pipe open."
tags: [Science, Mathematics, Programming, subprocess, minimisation]
comments: false
math: true
---

## The Problem:

We want to programmatically find the minimum value of the following simple function

\\[
 f(x,y) = x^2 + y^2 + 3
 \tag{1}
\\]
We already know the answer to this, \\(\text{min}(f(x,y)) = 3 \\) when \\(x=0,y=0\\) but the catch and purpose of this post is that the function declaration is written in C++ whilst we want to optimise it from Python. Why even bother with this? Programming the minimisation of such a simple function is trivial in both languages but there are several reasons why this process might be very useful for more complex problems, including:

 1. Limitations of the Python language to solve a problem, requiring lower level languages.

 2. You want to use Python as a front-end. Python is excellent at 'gluing' other languages together, meaning you optimise the computationally demanding code in lower level languages such as C++ or Fortran with everything else done in Python saving a lot of programming time.

 3. The quantity being minimised takes too long to calculate in Python and is much faster in C++. 

Conventionally, calling C++ from Python is achieved by wrapping the C++ program so as to directly call it from Python. There are multiple ways to achieve this including  [Boost Python](https://www.boost.org/doc/libs/1_66_0/libs/python/doc/html/index.html), [Swig](http://www.swig.org/) and [Cython](https://cython.org/) all of which work very well; but can require considerable time to accomplish depending on how complicated the program you are trying to wrap is. Time is precious in scientific research so rather than spend time wrapping a complicated C++ program I will sometimes use an alternative, **subprocessing**.

## Subprocessing
The [Subprocess](https://docs.python.org/2/library/subprocess.html) Python module is a way to spawn new processes and connect to their input/output pipes. This means that we can directly call the C++ executable from Python, pass it input and then grab the return code/returned result to Python. The benefit of this approach is that it is incredibly easy to implement. 


Let's think about the individual tasks C++ and Python need to accomplish to make them "talk" to one another in order to minimise the function in equation (1).

**C++:**
 1. Calculate the value of the function for specified values of \\(x\\) and \\(y\\).
 
 2. Output this value so Python is able to grab it.
 
 3. Wait for Python to provide new values of \\(x\\) and \\(y\\). 
 
 4. Repeat until the Python program terminates the executable.

**Python:**
 
 1. Control the minimisation process, easily done using `scipy`.
 
 2. Initiate the subprocess to call the C++ executable.
 
 3. Pass the values of \\(x\\) and \\(y\\) for each optimisation step to the C++ subprocess.
 
 4. Repeat until convergence tolerance is achieved.

 5. Kill the subprocess so it does not hang indefinitely.

## C++ Program

The C++ program, `Subprocess.cpp` looks like this

```cpp
#include <iostream>
#include <math.h>

// Define the function to be minimised, x**2 + y**2 + 3
double FuncValue (double x, double y)
{    
    return pow(x,2) + pow(y,2) + 3;
}

int main()
{	
    // Declare the data types of the input values x and y
    double x;
    double y;

    std::printf("This C++ program will be called from Python continuously\n");
    do 
    {   
        // This reads in the command line input and assigns to x and y
        std::cin >> x; 
        std::cin >> y; 
        
        // Print function value to the terminal
        std::printf("%.16f \n", FuncValue(x, y));
        
        // Clear cin
        std::cin.clear();
    } while (1 < 2); // Keep do loop running indefinitely

    return 0;
}
```
Note, the newline character `\n` is **vital** in the `printf` statements. This sends the signal to Python that something has been returned. if they are not included the program will hang and nothing will be returned. For this example the Python program will terminate the executable by using its session id, which is not very elegant but works very well. A more elegant option is to specify a custom flag for the `while` condition. For example `while (std::cin.get() != 'KILL')`. This then means that we can pass the word `KILL` from Python and the executable terminates itself. 

An important feature of this program is the `do while` loop. Inside this loop `cin` is declared which means the program will continually loop waiting for user input, which will be given to it by Python. We have to make it wait otherwise we would need to terminate and restart the executable for every optimisation step. This is no issue for such a trivial example, but adds substantial computational overhead for more complex examples. 

To compile this program type the following into the terminal

```bash
g++ Subprocess.cpp -o Subprocess
```
## Python Program
The Python program, `Subprocess.py` looks like this

```python
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

    print("x={:.6E} y={:.6E} f(x,y)={:.16f}".format(x, y, FuncValue))

    return FuncValue

# Initiate the subprocess by calling the Subprocess executable we just compiled
# 'universal_newlines' is a handy option as it means we don't have to supply the input as a byte string or decode the output string as utf-8 as the file objects stdout and stderr are opened as text files in universal newlines mode
# preexec_fn=os.setsid assigns a session id to the subprocess which we use later on to terminate it
process = Popen('./Subprocess', stdin=PIPE, stdout=PIPE, universal_newlines=True, shell=True, preexec_fn=os.setsid)

# Tell Python that the first line output by C++ is: "This C++ program will be called from Python continuously" otherwise it will prevent the values of x and y being piped
print(process.stdout.readline()) 

# provide the starting values for x and y
initial_guess = [2,2]

# Minimise the function CallCpp using the Nelder-Mead algorithm
Func_min = minimize(CallCpp, initial_guess, method='nelder-mead')

# Kill the subprocess
os.killpg(os.getpgid(process.pid), signal.SIGTERM)
```
The newline character in `process.stdin.write('{} {}\n'.format(x, y))` is also **vital**. If this is not included then we have not actually entered the \\(x\\) and \\(y\\) values into the executable. Think of the `\n` acting as the enter key.
## The Result
To run the Python program type the following
```bash
python3 Subprocess.py
```
This gives the following output
```text
This C++ program will be called from Python continuously

x=2.000000E+00 y=2.000000E+00 f(x,y)=11.0000000000000000

x=2.100000E+00 y=2.000000E+00 f(x,y)=11.4100000000000001

x=2.000000E+00 y=2.100000E+00 f(x,y)=11.4100000000000001

x=2.100000E+00 y=1.900000E+00 f(x,y)=11.0199999999999978

x=2.000000E+00 y=1.900000E+00 f(x,y)=10.6099999999999977

x=1.950000E+00 y=1.850000E+00 f(x,y)=10.2249999999999961

x=1.850000E+00 y=1.950000E+00 f(x,y)=10.2249999999999979

x=1.800000E+00 y=1.800000E+00 f(x,y)=9.4799999999999933

x=1.700000E+00 y=1.700000E+00 f(x,y)=8.7799999999999923

x=1.800000E+00 y=1.600000E+00 f(x,y)=8.7999999999999901

x=1.550000E+00 y=1.450000E+00 f(x,y)=7.5049999999999875

x=1.350000E+00 y=1.250000E+00 f(x,y)=6.3849999999999847

x=1.250000E+00 y=1.350000E+00 f(x,y)=6.3849999999999865

x=9.000000E-01 y=9.000000E-01 f(x,y)=4.6199999999999850

x=5.000000E-01 y=5.000000E-01 f(x,y)=3.4999999999999885

x=6.000000E-01 y=4.000000E-01 f(x,y)=3.5199999999999871

x=-2.500000E-01 y=-3.500000E-01 f(x,y)=3.1850000000000098

x=-1.050000E+00 y=-1.150000E+00 f(x,y)=5.4250000000000522

x=-3.500000E-01 y=-2.500000E-01 f(x,y)=3.1850000000000112

x=-1.100000E+00 y=-1.100000E+00 f(x,y)=5.4200000000000523

x=1.000000E-01 y=1.000000E-01 f(x,y)=3.0199999999999969

x=2.000000E-01 y=-4.551914E-15 f(x,y)=3.0399999999999956

x=5.500000E-01 y=4.500000E-01 f(x,y)=3.5049999999999879

x=-5.000000E-02 y=-1.500000E-01 f(x,y)=3.0250000000000026

x=-1.500000E-01 y=-5.000000E-02 f(x,y)=3.0250000000000039

x=-6.250000E-02 y=-3.750000E-02 f(x,y)=3.0053125000000018

x=8.750000E-02 y=2.125000E-01 f(x,y)=3.0528124999999964

x=-1.562500E-02 y=-5.937500E-02 f(x,y)=3.0037695312500010

x=-1.781250E-01 y=-1.968750E-01 f(x,y)=3.0704882812500061

x=3.046875E-02 y=2.578125E-02 f(x,y)=3.0015930175781240

x=7.734375E-02 y=3.906250E-03 f(x,y)=3.0059973144531233

x=-2.753906E-02 y=-2.714844E-02 f(x,y)=3.0014954376220713

x=1.855469E-02 y=5.800781E-02 f(x,y)=3.0037091827392568

x=1.000977E-02 y=2.866211E-02 f(x,y)=3.0009217119216913

x=-4.799805E-02 y=-2.426758E-02 f(x,y)=3.0028927278518691

x=1.085205E-02 y=1.326904E-02 f(x,y)=3.0002938345074650

x=4.840088E-02 y=6.907959E-02 f(x,y)=3.0071146348118765

x=-8.554077E-03 y=-3.091431E-03 f(x,y)=3.0000827291794123

x=-7.711792E-03 y=-1.848450E-02 f(x,y)=3.0004011483676734

x=-3.281403E-03 y=-6.697845E-03 f(x,y)=3.0000556287367366

x=-2.268753E-02 y=-2.305832E-02 f(x,y)=3.0010464101203258

x=2.467155E-03 y=4.187202E-03 f(x,y)=3.0000236195204342

x=7.739830E-03 y=5.807877E-04 f(x,y)=3.0000602422829976

x=3.666353E-03 y=-3.372669E-04 f(x,y)=3.0000135558949523

x=9.414911E-03 y=1.054778E-02 f(x,y)=3.0001998962380498

x=-1.073241E-04 y=-2.386439E-03 f(x,y)=3.0000057066088361

x=1.091874E-03 y=-6.910908E-03 f(x,y)=3.0000489528405141

x=2.123335E-03 y=1.412675E-03 f(x,y)=3.0000065042015858

x=-1.650342E-03 y=-6.364971E-04 f(x,y)=3.0000031287584656

x=-4.308690E-03 y=-7.861122E-04 f(x,y)=3.0000191827831104

x=-3.881001E-03 y=-4.435611E-03 f(x,y)=3.0000347368153428

x=6.222509E-04 y=-4.939660E-05 f(x,y)=3.0000003896361882

x=-9.207673E-04 y=1.700545E-03 f(x,y)=3.0000037396661341

x=-7.174065E-04 y=6.787991E-04 f(x,y)=3.0000009754403734

x=1.555187E-03 y=1.265900E-03 f(x,y)=3.0000040211075860

x=-8.489601E-04 y=-1.608979E-04 f(x,y)=3.0000007466213736

x=4.906973E-04 y=-8.890937E-04 f(x,y)=3.0000010312714123

x=-4.153806E-04 y=2.868259E-04 f(x,y)=3.0000002548101241

x=1.055830E-03 y=3.983273E-04 f(x,y)=3.0000012734424444

x=-3.727625E-04 y=-2.109164E-05 f(x,y)=3.0000001393967111

x=-1.410394E-03 y=3.151309E-04 f(x,y)=3.0000020885184679

x=1.140897E-04 y=4.173527E-05 f(x,y)=3.0000000147582888

x=1.567078E-04 y=-2.661823E-04 f(x,y)=3.0000000954103454

x=6.435599E-04 y=-2.033554E-04 f(x,y)=3.0000004555228053

x=-1.186819E-04 y=-6.665758E-05 f(x,y)=3.0000000185286173

x=-1.613000E-04 y=2.412600E-04 f(x,y)=3.0000000842240597

x=-8.179803E-05 y=1.143994E-04 f(x,y)=3.0000000197781436

x=7.720585E-05 y=-1.393217E-04 f(x,y)=3.0000000253712864

x=-4.204706E-05 y=5.096913E-05 f(x,y)=3.0000000043658073

x=1.907245E-04 y=1.593620E-04 f(x,y)=3.0000000617720701

x=-4.133027E-05 y=-1.015269E-05 f(x,y)=3.0000000018112689

x=-1.974670E-04 y=-9.188273E-07 f(x,y)=3.0000000389940675

x=3.620051E-05 y=3.107174E-05 f(x,y)=3.0000000022759301

x=3.691729E-05 y=-3.005008E-05 f(x,y)=3.0000000022658937

x=-4.061349E-05 y=-7.127451E-05 f(x,y)=3.0000000067295112

x=1.699701E-05 y=5.485181E-06 f(x,y)=3.0000000003189857

x=-6.125056E-05 y=2.538257E-05 f(x,y)=3.0000000043959054

x=1.237533E-05 y=-1.619192E-05 f(x,y)=3.0000000004153269

x=7.070261E-05 y=-5.540448E-07 f(x,y)=3.0000000049991664

x=-1.332205E-05 y=-7.753028E-06 f(x,y)=3.0000000002375864

x=-8.700374E-06 y=1.392407E-05 f(x,y)=3.0000000002695764

Optimization terminated successfully.
         Current function value: 3.000000
         Iterations: 45
         Function evaluations: 81

```
Resulting in a minimum function value of \\(f(x,y) = 3\\) where \\(x = -8.700374\times 10^{-6}\\) and \\(y=1.392407 \times 10^{-5}\\), which are numerically close to 0 given the convergence tolerance specified.
# Conclusion

We have minimised equation (1) using a combination of Python and C++. Python controlled the minimisation by piping the parameters, \\(x\\) and \\(y\\) to the C++ executable which then returned the function value for Python to decide the values of \\(x\\) and \\(y\\) for the next minimisation step. This continued until the convergence tolerance was achieved at which point Python terminated the subprocess and both programs ended.

It is advised to wrap a C++ program using conventional methods, but sometimes it is far easier to use a subprocess achieving the same results.