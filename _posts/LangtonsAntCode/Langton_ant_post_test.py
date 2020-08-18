import numpy as np
import itertools
from matplotlib import pyplot as plt
from matplotlib import animation, colors

square_colours = ['white', 
                  'black', 
                  'indianred', 
                  'dodgerblue'
                 ]
# Set the dimensions of the grid
dim = 200
# Set the number of steps the ant should take
ant_steps = 1000000
# Build a corresponding numpy array of dimensions (dim,dim)
grid = np.array(np.zeros((dim, dim)))
# Define a variable to represent the current ant_position of the ant on the board
ant_pos = np.array([[int(dim/2)], [int(dim/2)]])
# Define a variable to represent the direction ant is currently moving
direction = np.matrix([[1], [0]])

# Clockwise rotation matrix
clockwise_rot     = np.array([[0, 1], [-1, 0]])
# Anti-clockwise rotation matrix
anticlockwise_rot = np.array([[0, -1], [1, 0]])

def move_ant(grid, ant_pos, direction):
    '''
    Controls the movement of the ant by a single square

    This function takes the current position and the direction of the ant and updates it via the 2 rules specified above as it takes its next stdep. It then updates the new position, direction and square colour for the next step.

    Parameters:
    grid (np.array)      : This is the grid of dimension, dim, that the ant moves around on
    ant_pos (np.array)   : This represents the ants' position defined as a numpy array of its x,y coordinate on the grid
    direction(np.matrix) : This represents the direction that the ant will move in on this step. 

    Returns:
    None: No explicit return 
    '''
    # Note ant_pos[0, 0], ant_pos[1, 0] is the x and y coordinate of the ant on the grid
    ant_pos[:] = ant_pos + direction
    x_ant_pos = ant_pos[0, 0]
    y_ant_pos = ant_pos[1, 0]

    # LRRL instance
    if any(i == dim or i == 0 for i in ant_pos):
        print("Hit the edge of the board!")
        exit()
    elif grid[x_ant_pos, y_ant_pos] == 0:            # Landed on white
        direction[:] = anticlockwise_rot * direction # Ant turns left and moves forward
        grid[x_ant_pos, y_ant_pos] = 1               # As landed on white square, change to black square
    elif grid[x_ant_pos, y_ant_pos] == 1:            # Landed on black
        direction[:] = clockwise_rot * direction     # Ant turns right and moves forward
        grid[x_ant_pos, y_ant_pos] = 2               # As landed on black square, change to red square
    elif grid[x_ant_pos, y_ant_pos] == 2:            # Landed on red
        direction[:] = clockwise_rot * direction     # Ant turns right and moves forward  
        grid[x_ant_pos, y_ant_pos] = 3               # As landed on red square, change to blue square
    else:                                            # Landed on blue 
        direction[:] = anticlockwise_rot * direction # Ant turns left and moves forward  
        grid[x_ant_pos, y_ant_pos] = 0               # As landed on blue square, change to white square

fig = plt.figure()
ax = fig.add_subplot(111)

# Define custom discrete colour map
cmap = colors.ListedColormap(square_colours)

boundaries = [0, 0, 1, 1, 2, 2, 3, 3]
            # |w |  |bl|  |r |  |b |

norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
# Define the 'image' to be created using imshow. Specify the custom colour map
im = plt.imshow(grid, cmap=cmap, norm=norm)

for i in range(ant_steps):
    move_ant(grid, ant_pos, direction) 


im.set_data(grid)

# Hide the axes
ax.axis('off')
# Report the number of ant steps taken in the title of the image
plt.suptitle("No Steps Taken = {:d}".format(ant_steps), fontsize=13)
# Show the image
plt.show()