import numpy as np
import itertools
from matplotlib import pyplot as plt
from matplotlib import animation, colors
import seaborn as sb

# Define the colour palette to be used
colours = "custom1"

if colours is "seaborn":
    square_colours = sb.color_palette(palette='deep') # Options: 'deep', 'muted', 'pastel', 'bright', 'dark', and 'colorblind'
elif colours is "custom1":
    square_colours = ['white', 'black', 'indianred', 'dodgerblue', 'gold', "steelblue", "tomato", "slategray", "plum", "seagreen", "gray"]
elif colours is "custom2":
    square_colours = ["#000000", "#FF0000", "#444444", "#FFFF00", "#00FF00", "#00FFFF",  "#0000FF", "#9900FF"]
else:
    print("Colour palette not recognised. Exiting")
    exit()

'''
User Changeable Input
''' 
# Set the dimensions of the grid
dim = 700
# Set the number of steps the ant should take
ant_steps = 1000000
# Tell the program what moveset to give the ant
ant_move_list = 'LRRRRLLLRRR'

'''
End User Changeable Input
'''

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

# Calculate the length of the ant_move_list
len_ant_move_list = len(ant_move_list)
# Extract the unique 'L' and 'R' moves in order
unique_moves = [i for i in ant_move_list]
# Assign the correct rotation matrix to each 'L' and 'R' choice and store in a list
rot_matrices = [anticlockwise_rot if i == 'L' else clockwise_rot for i in unique_moves]
# Assign each of these unique 'L' and 'R' letters to a discrete integer which will later represent a colour 
colour_indices = [i for i in range(len(ant_move_list))]

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
    # Create the next ant position
    ant_pos[:] = ant_pos + direction
    # Extract the x coordinate of this new position
    x_ant_pos = ant_pos[0, 0]
    # Extract the y coordinate of this new position
    y_ant_pos = ant_pos[1, 0]

    if any(i == dim or i == 0 for i in ant_pos):
        print("Hit the edge of the board!")
        exit()
    if grid[x_ant_pos, y_ant_pos] in colour_indices:
        index = grid[x_ant_pos, y_ant_pos]
        grid[x_ant_pos, y_ant_pos] = colour_indices[int(index+1) % len(colour_indices)]
        direction[:] = rot_matrices[int(index)] * direction
    else:
        print("Index not in colour indices. Exiting")
        exit()

# Define the boundaries for each discrete colour
# The form is [0, 0, 1, 1, 2, 2, 3, 3]
#               c1    c2    c3    c4
# where cX represents a unique colour 
boundaries = list(itertools.chain.from_iterable(itertools.repeat(x, 2) for x in colour_indices))

fig = plt.figure()
ax = fig.add_subplot(111)

# Define custom discrete colour map
cmap = colors.ListedColormap(square_colours[0:len_ant_move_list])

norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
# Define the 'image' to be created using imshow. Specify the custom colour map
im = plt.imshow(grid, cmap=cmap, norm=norm)

for i in range(ant_steps):
    move_ant(grid, ant_pos, direction) 

im.set_data(grid)
# Hide the axes
ax.axis('off')
# Report the movement pattern number of ant steps taken in the title of the image
plt.suptitle("{}, No. Steps Taken = {:d}".format(ant_move_list, ant_steps), fontsize=13)
# Save the image

plt.savefig("{}_{}.png".format(ant_move_list, ant_steps), format='png')
# Show the image
#plt.show()