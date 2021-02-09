import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import matplotlib.pyplot as plt

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
    '''
    This function is iterated over by FuncAnimation to produce the frames necessary to produce the animation 
    
    Parameters:
    -----------
    i : int
        This is the index used to iterate through the coordinate arrays
    
    Returns:
    --------
    line : mpl_toolkits.mplot3d.art3d.Line3D
           This contains the line information which is updated each time the function is iterated over
    '''
    line.set_data(X[:i], Y[:i])
    line.set_3d_properties(Z[:i])

    return line

anim = FuncAnimation(fig, animate, interval=10, frames=N, repeat=False)

HTML(anim.to_html5_video()) # This is just for the Jupyter notebook which will not show the animation when using plt.show()

# if running locally, uncomment these two lines and remove the HTML command above.
anim.save('BrownianMotion.mp4', writer='ffmpeg')