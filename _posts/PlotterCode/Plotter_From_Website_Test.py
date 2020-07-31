#!/usr/bin/env python3
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import PySimpleGUI as sg
import sys
import os
import itertools
# Set the seaborn dark grid styling
sns.set_style("darkgrid")
# Set the theme for PysimpleGUI
sg.theme('DarkAmber')   

if len(sys.argv) == 1:
    event, fnames = sg.Window('Select File(s) you wish to plot.').Layout([[sg.Text('Note, select multiple files by holding ctrl and clicking the number required.')],
                                                                          [sg.Input(key='_FILES_'), sg.FilesBrowse()], 
                                                                          [sg.OK(), sg.Cancel()]]).Read(close=True)
    # Close the window if cancel is pressed
    if event in (sg.WIN_CLOSED, 'Cancel'):
        exit()

else:
    # Check to see if any files were provided on the command line
    fnames = sys.argv[1]

# If no file names are selected, exit the program
if not fnames:
    sg.popup("Cancel", "No filename supplied")
    raise SystemExit("Cancelling: no filename supplied")

# Count the number of files provided
fnames = fnames['_FILES_'].split(';')
no_files = len(fnames)

# List the available colours for the plots
matplotlib_colours = ["dodgerblue", "indianred", "gold", "steelblue", "tomato", "slategray", "plum", "seagreen", "gray"]
# List the line-styles you want
matplotlib_linestyles = ["solid", "dashed", "dashdot", "dotted"]

headings = ['X,Y INDICES', '  TYPE', 'COLOUR','LINE', '  LEGEND']  # the text of the headings

# Create the layout of the Window
layout = [  [sg.Text('You can use LaTeX math code for axis labels and legend entries, e.g. $\mathbf{r}$', font=('Courier', 10))],
            [sg.Text('To use regular text in math mode use $\mathrm{Text}$\n')],
            [sg.Text('_'  * 100, size=(100, 1))], # Add horizontal spacer  
            [sg.Text('X-axis label:'), 
              sg.InputText('')],
            [sg.Text('Y-axis label:'), 
              sg.InputText('')],
            [sg.Text('_'  * 100, size=(100, 1))], # Add horizontal spacer  
            [sg.Text('                                        ')] + [sg.Text(h, size=(11,1)) for h in headings],  # build header layout
            *[[sg.Text('File: {}'.format(os.path.basename(os.path.normpath(i))), size=(40, 1)), 
               sg.InputText('X', size=(5, 1)),
               sg.InputText('Y', size=(5, 1)),
               sg.InputCombo(values=('point', 'line')),
               sg.InputCombo(values=(matplotlib_colours)),
               sg.InputCombo(values=(matplotlib_linestyles)),
               sg.InputText('Enter Legend Label', size=(20, 1)),
              ] for i in fnames
             ],
            [sg.Text('_'  * 100, size=(100, 1))], # Add horizontal spacer
            [sg.Button('Plot'), sg.Button('Cancel')],
         ]

# Create the Window
window = sg.Window('Plot v1-01', layout)
# Read in the events and values       
event, values = window.read()
# If cancel is pressed then close the window and exit
if event in (sg.WIN_CLOSED, 'Cancel'):
    exit()

window.close()

# Access the values which were entered and store in lists
xAxisLabel = values[0]
yAxisLabel = values[1]

listOfDataSets = []
legendLabels   = []
xcols          = []
ycols          = []
cols_to_use    = []
plot_type      = []
plot_colour    = []
plot_line      = []
i = 2
for file in fnames:
    # Append the data files to a single list
    listOfDataSets.append(file)
    # Append the column indices to a list for later
    xcolindex = int(values[i]) - 1 # index 2
    i += 1
    ycolindex = int(values[i]) - 1 # index 3
    # Append the separate x and y column indices to their respective lists. These are used when plotting using Seaborn below
    xcols.append(xcolindex)
    ycols.append(ycolindex)
    # Append both the x and y to a combined list in order to construct the DataFrame object
    cols_to_use.append([xcolindex, ycolindex])
    # Append the type of plot [ scatter | line ]
    i += 1
    plot_type.append(values[i]) # index 4
    # Append the colour of the plot
    i += 1
    plot_colour.append(values[i]) # index 5
    # Append the linestyle of the plot
    i += 1
    plot_line.append(values[i]) # index 6
    # Append the user specified legend labels to a list for later
    i += 1
    legendLabels.append(values[i]) # index 7
    i += 1

fig, ax = plt.subplots(figsize=(4, 4))

# Interate over the colours and line styles provided from the GUI 
plot_colour = itertools.cycle(plot_colour)
plot_line = itertools.cycle(plot_line)

for i in range(0, len(listOfDataSets)):
    # Read in data set into a pandas dataframe. Note [cols_to_use[i]] at the end maintains column index order 
    df = pd.read_csv(listOfDataSets[i], usecols=cols_to_use[i], sep="\s+|\t+|\s+\t+|\t+\s+|,\s+|\s+,", header=None, engine='python')[cols_to_use[i]]
    # Plot the data set using Seaborn and set legend labels from user specified ones above
    if plot_type[i] == 'point':
        ax.scatter(df[xcols[i]], df[ycols[i]], color=next(plot_colour), s=10, label=r'{}'.format(legendLabels[i]))
    elif plot_type[i] == 'line':
        ax.plot(df[xcols[i]], df[ycols[i]], color=next(plot_colour), linestyle=next(plot_line), label=r'{}'.format(legendLabels[i]))
    else:
        # If a plot type is not specified [point | line] then default to line
        print("\n No option for plot type specified, defaulting to line plot")
        ax.plot(df[xcols[i]], df[ycols[i]], color=next(plot_colour), linestyle=next(plot_line),  label=r'{}'.format(legendLabels[i]))

    # Work out the minimum and maximum values in the columns to get the plotting range correct
    xmin = df[xcols[i]].min()
    xmax = df[xcols[i]].max()
    ymin = df[ycols[i]].min()
    ymax = df[ycols[i]].max()
    # Set axis limits
    plt.xlim(xmin, None)
    plt.ylim(ymin, None)
    # Set the x and y axis labels from the user specified ones above
    plt.xlabel(r'{}'.format(xAxisLabel))
    plt.ylabel(r'{}'.format(yAxisLabel))
    # Show the legend
    plt.legend()

# Finally show the plot on screen
plt.show() 