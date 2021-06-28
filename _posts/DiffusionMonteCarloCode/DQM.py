import numpy as np
import math
import matplotlib.pyplot as plt
import copy
from mpl_toolkits.mplot3d import Axes3D

numAtoms = 1
initialWalkers = 1000
wvnmbr = 4.55634e-6
omega = 2000.0000*wvnmbr #in atomic units
dimensions = 1
mass = 1836.35/2

# imaginary time step
deltaT = 5.0 #simulation parameter - adjustable

alpha = 1.0/(2.0*deltaT) #simulation parameter - adjustable

D = 0.5
sigma = math.sqrt((2*D*deltaT)/mass)

class Walker:
    def __init__(self):
        self.coords = np.zeros((numAtoms,dimensions)) #1d surface
        self.weight = 0.0
        self.WalkerV = 0.0

myWalkers = [Walker() for r in range(initialWalkers)]

def getPotentialForWalkers(): #use coordinates of walkers to get V
    omsqd = math.pow(omega,2)
    prefactor = 0.50000*mass*omsqd
    for d in range(len(myWalkers)):
        crds = float(myWalkers[d].coords)
        crdssq = math.pow(crds,2)
        myWalkers[d].WalkerV = (prefactor*crdssq)

def giveBirth(singleWalker):
    myWalkers.append(singleWalker)

def deleteWalker(walkerIndex):
    del myWalkers[walkerIndex]

def birthOrDeath(vref):
    deathAr = []
    birthAr = []
    for y in range(len(myWalkers)):
        Rng = np.random.random()
        curV = myWalkers[y].WalkerV
        if curV > Vref:
            P = -1*(curV - Vref)*deltaT
            exP = math.exp(P)
            if exP < Rng:
                deathAr.append(y)
        else:
            P = -1*(curV - Vref) * deltaT
            exP = math.exp(P) - 1
            if exP > Rng:
                singleWalker = copy.deepcopy(myWalkers[y])
                birthAr.append(singleWalker)

    if deathAr: #If it's not empty
        for k in reversed(deathAr):
            deleteWalker(k)
    if birthAr:
        for j in birthAr:
            giveBirth(j)
    del deathAr[:]
    del birthAr[:]

def moveRandomly():
    #choose a random number from gaussian distribution (1/2pisig)(e^(-dx^2/1sig^2))
    for p in range(len(myWalkers)):
        gaussStep = np.random.normal(loc = 0.0, scale=sigma)
        myWalkers[p].coords = myWalkers[p].coords + gaussStep

def getVref(): #Use potential of all walkers to calculate vref
    varray = np.array([k.WalkerV for k in myWalkers])
    Vbar = np.average(varray)
    alphaterm = alpha*((len(myWalkers)-initialWalkers)/initialWalkers)
    vref = Vbar - alpha*((float(len(myWalkers))-float(initialWalkers))/float(initialWalkers))
    print("AlphaTerm ! = ",alphaterm)
    print("Numwalkers", len(myWalkers))
    return vref

#Start!
vrefAr = np.zeros(1000)
xc = []
popAr = np.zeros(1000)
Vref = 10000
print("initialVref = ",Vref)

for i in range(1000):
    moveRandomly()
    getPotentialForWalkers()

    print("Beginning VRef = ",Vref)
    print("Beginning NWalkers",len(myWalkers))
    if i==0:
        Vref = getVref()
    # if i>=975:
    #     crds = [float(o.coords) for o in myWalkers]
    #     pots = [c.WalkerV for c in myWalkers]
    #     plt.scatter(crds,pots)
    #     plt.plot(np.linspace(-2,2,100),[Vref]*100)
    #     plt.show()
    birthOrDeath(Vref)
    #After the purge/birth cycle, Get the new potential

    getPotentialForWalkers()
    #This will allow us to get a new vref for next cycle.
    Vref = getVref()
    print("Average potential (AU)= ", np.average([p.WalkerV for p in myWalkers]))
    print("Post Birth Vref", Vref)
    vrefAr[i] = Vref
    print("Post Birth NWalkers = ",len(myWalkers))
    print("End Loop", i + 1)
    popAr[i] = len(myWalkers)
    xc.append([float(n.coords) for n in myWalkers])
    print("Average coordinates = ", np.average([float(n.coords) for n in myWalkers]))

fff,ax1 = plt.subplots()
ax1 = plt.axes(projection='3d')
b=1
for q in xc:
    r,binz = np.histogram(q,bins=100,range=[-10,10])
    ax1.plot(np.linspace(-10,10,100),r,b)
    b+=1
plt.show()

var = [h for h in vrefAr]
x = np.linspace(0,1000,1000)
plt.plot(x,vrefAr)
om = [omega/2]*len(x)
plt.plot(x,om)
plt.show()
plt.plot(x,popAr)
plt.show()