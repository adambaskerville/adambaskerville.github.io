import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

def energyplot(x):
    '''
    The hydrogenic energy expression
    '''
    return (x**2)/2 - x

def energyplot_gaussian(x):
    '''
    The hydrogenic energy expression
    '''
    #return x - x**0.5
    return (3/2)*x - (2*(2**0.5)/(np.pi**0.5))*(x**0.5)

x = np.linspace(0, 3, 100)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

#spine placement data centered
ax.spines['left'].set_position(('data', 0.0))
ax.spines['bottom'].set_position(('data', 0.0))
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
# Set position of x and y axis labels
ax.xaxis.set_label_coords(1, 0.33)
ax.yaxis.set_label_coords(0.08, 1)

sns.set_context("talk", font_scale=1.5, rc={"lines.linewidth": 2.5})

plt.xticks(np.arange(0, 3, 1.0))
plt.yticks(np.arange(-0.5, 2, 1.0))
# Custom x axis label
plt.xlabel(r"$c$", fontsize=15)
# Custom y axis label
ylab = plt.ylabel(r"$E$", fontsize=15)
# Rotate the y axis label
ylab.set_rotation(0)

# Plot the function
plt.plot(x, energyplot(x), label=r"$E(c)_{\psi_1}$")
plt.plot(x, energyplot_gaussian(x), label=r"$E(c)_{\psi_2}$")
plt.legend(fontsize=12, loc = 'upper center')
# Show the plot
plt.show()