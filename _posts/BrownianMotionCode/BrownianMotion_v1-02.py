import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import mpl_toolkits.mplot3d.axes3d as p3
import seaborn as sns
sns.set_style("darkgrid", {'axes.grid' : False})

## Note look up Wiener process when writing blog post

np.random.seed(5)

fig = plt.figure()
ax = p3.Axes3D(fig)

N = 1000 # Number of points
T = 1.0
dt = T/(N-1)

dX = np.sqrt(dt) * np.random.randn(1, N)
X = np.cumsum(dX) # Cumulatively add these values to get the x positions of the particle

dY = np.sqrt(dt) * np.random.randn(1, N) 
Y = np.cumsum(dY) # Cumulatively add these values to get the y positions of the particle

dZ = np.sqrt(dt) * np.random.randn(1, N) 
Z = np.cumsum(dZ) # Cumulatively add these values to get the z positions of the particle

line, = ax.plot(X, Y, Z, color='blue', lw=1)

def animate(i):
    line.set_data(X[:i], Y[:i])
    line.set_3d_properties(Z[:i])

    return line

anim = FuncAnimation(fig, animate, interval=10, frames=N, repeat=False)

plt.draw()

# Plot start and end points
#ax.scatter(X[0], Y[0], Z[0], 'ro')
#ax.scatter(X[-1], Y[-1], Z[-1], 'yo')

plt.show()