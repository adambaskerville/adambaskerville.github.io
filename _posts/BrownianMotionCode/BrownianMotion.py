import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import seaborn as sb

np.random.seed(5)

fig = plt.figure()
ax = plt.axes(xlim=(-2, 2), ylim=(-2, 2))
line, = ax.plot([], [], lw=2, color='blue')

N = 1000 # Number of points
T = 1.0
dt = T/(N-1)

dX = np.sqrt(dt) * np.random.randn(1, N)
X = np.cumsum(dX) # Cumulatively add these values to get the x positions of the particle

dY = np.sqrt(dt) * np.random.randn(1, N) 
Y = np.cumsum(dY) # Cumulatively add these values to get the y positions of the particle

dZ = np.sqrt(dt) * np.random.randn(1, N) 
Z = np.cumsum(dZ) # Cumulatively add these values to get the y positions of the particle

line = ax.plot(X, Y, color='blue', lw=2)[0]

def init():
    line.set_data([], [])
    return line,

def animate(i):
    line.set_data(X[:i], Y[:i])

    return line

anim = FuncAnimation(fig, animate, init_func=init, interval=10, frames=N, blit=False)
 
plt.draw()
plt.show()