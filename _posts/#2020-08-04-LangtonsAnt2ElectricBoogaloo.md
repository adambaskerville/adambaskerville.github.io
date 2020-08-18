---
layout: post
title: "T>T: Langton's Ant 2: Electric Boogaloo"
date: 2020-07-26
excerpt: "Further investigation into the complex emergent behaviour of langton's ant."
tags: [science, mathematics, programming, matplotlib, script, langtons, ant, python]
comments: false
math: true
---

In a [previous post](https://adambaskerville.github.io/posts/LangtonsAnt/) we programmed a simple implementation of Langton's ant and before reading ahead, it is recommended to first read this initial post. Our ant obeyed the following movement pattern:

1) If you are on a white square, turn 90° clockwise, flip the colour of the square, then move forward one square.

2) If you are on a black square, turn 90° anti-clockwise, flip the colour of the square, then move forward one square.

This movement pattern is classed as the `RL` scheme, as the first step involves turning right, whilst the second step involves turning left, alternating the colour each time. There are an infinity of possible extensions to this idea; and in this post we will explore the instance of multiple colours from different step patterns. 

By the end we will have a generalised program for any movement pattern, e.g. `LRRL`

# `LRRL` Langton's Ant

A good tip before attempting to create a completely generalised program is to program several individual test cases so you fully understand the required behaviour; looking for any pitfalls or edge cases which need special consideration. The previous post handled the `RL` movement scheme so let's consider a slightly more complicated one, the `LRRL` scheme described as follows:

1. If you are on a white square, turn 90° anti-clockwise, flip the colour of the square to <span style="color:black">black</span>, then move forward one square. (`L`)

2. If you are on a <span style="color:black">black</span> square, turn 90° clockwise, flip the colour of the square to <span style="color:red">red</span>, then move forward one square. (`R`)

3. If you are on a <span style="color:red">red</span> square, turn 90° clockwise, flip the colour of the square to <span style="color:blue">blue</span>, then move forward one square. (`R`)

4. If you are on a <span style="color:blue">blue</span> square, turn 90° anti-clockwise, flip the colour of the square to white, then move forward one square. (`L`)

This looks complicated but the ant is still just turning left and right, only this time it has two more colours to paint the squares. Consider the pseudo-code of the `LRRL` scheme:

```
if ant on white square:
    turn left, change square to black, move forward one step
elif ant on black square:
    turn right, change square to red, move forward one step
elif ant on red square:
    turn right, change square to blue, move forward one step 
elif ant on blue square:
    turn left, change square to white, move forward one step 
```
In the previous post we used `0` to represent white and `1` to represent black so we will add more colours by just incrementing this index; red becomes `2` and blue becomes `3`. We can store these colours in a list for easy access later on:

```python
square_colours = ['white', 
                  'black', 
                  'indianred', 
                  'dodgerblue'
                 ]
```
We will now convert the above pseudo-code into Python code using the same syntax from the first Langton's ant post seen in the `move_ant` function. Read the comments for a description of each conditional statement:

```python
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
```
We now have our colours and movement pattern, so the next step is to define the colour map used by `matplotlib`. We need a discrete colour scheme with well defined boundaries for each colour so there is no gradient across colours as each square is **one** colour. First, assign the list of square colours to a colour map object using `ListedColormap` defining the colours which the colour map will be built from:

```python
cmap = colors.ListedColormap(square_colours)
```
Next we need to program the boundaries for each colour which correspond to the index of the colour in the list. We want the following

* White if grid value is 0
* <span style="color:black">Black</span> if grid value is 1
* <span style="color:red">Red</span> if grid value is 2
* <span style="color:blue">Blue</span> if grid value is 3

We can do this by setting a list of boundary values for each colour as follows:

```python
# Boundaries syntax:
# [lower_bound_colour1, upper_bound_colour1, lower_bound_colour2, upper_bound_colour2, ...]

boundaries = [0, 0, 1, 1, 2, 2, 3, 3]
            # |w |  |bl|  |r |  |b |
```
Each colour is now discrete, described by a single numerical value, so we can now set these boundaries using `BoundaryNorm` which will map the colour map indices to integer values.

```python
norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
```
We will copy over some of the code from the first post to start building our `LRRL` program, inserting our new set of rules inside the `move_ant` function. We are not going to copy over the `animate_ant` function as the animations can take a long time to generate and their file sizes can grow quickly. Instead we will just print the end result:

```python
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
dim = 50
# Set the number of steps the ant should take
ant_steps = 10000
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

# Create a figure object
fig = plt.figure()
ax = fig.add_subplot(111)

# Define custom discrete colour map
cmap = colors.ListedColormap(square_colours)
# Define the colour boundaries described above
boundaries = [0, 0, 1, 1, 2, 2, 3, 3]
            # |w |  |bl|  |r |  |b |
# Assign the boundaries using BoundaryNorm 
norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
# Define the 'image' to be created using imshow. Specify the custom colour map and the boundaries
im = plt.imshow(grid, cmap=cmap, norm=norm)

# Loop over the move_ant function for the number of steps required
for i in range(ant_steps):
    move_ant(grid, ant_pos, direction) 
# Set the final data for the image. Before we updated this each step for the animation. This time we only do it once at the end to print the final result
im.set_data(grid)

# Hide the axes
ax.axis('off')
# Report the number of ant steps taken in the title of the image
plt.suptitle("No Steps Taken = {:d}".format(ant_steps), fontsize=13)
# Show the image
plt.show()
```
This is a python implementation of the `LRRL` movement scheme for our ant. Here is a summary of what the program does:

1. We specified our square colours in the list `square_colours`.
2. The user specifies the number of steps the ant is to take along with the dimensions of the grid.
3. We specify a numpy array of zeros and place our ant at the mid point of the board.
4. We specify clockwise and anti-clockwise rotation matrices so our ant can turn right and left.
5. The main function `move_ant` controls everything for a single movement the ant makes. It reads the colour of the square as an integer value and decides what the next move should be, rotates the ant, changes the square colour and moves forward one square.
6. This is then looped over for the number of ant steps specified by the user.
7. The final numpy array of integers is printed to screen with each integer value set to its designated colour. The total number of steps is printed above it for reference. 

Now lets run the program to see if it works, and more importantly to see what path our ant takes!

