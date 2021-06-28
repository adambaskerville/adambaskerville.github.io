# GeneralDMC
General Diffusion Monte Carlo Package for solving the vibrational Schrodinger Equation.  Taking advantage of numpy vectorization, this is an efficient implementation of the algorithm described by James Anderson (https://doi.org/10.1063/1.432868).

This code is written in Python3, and it requires a potential energy surface to run for the system of interest.  This code can
uses all the cores available, but does not do multi-node calculations on HPCs.
