#Plots the orbits in real-time while a run is happening.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import itertools
import sys

file_name=str(sys.argv[1])
fos = open(file_name, 'r')
data = np.loadtxt(fos, delimiter=',')
x = data[:,0]   #time
y = data[:,1]   #energy

fig, ax = plt.subplots()
line, = ax.plot([], [], 'k-')
ax.margins(0.05)
ax.set_yscale('log')

def init():
    line.set_data(x[:2],y[:2])
    return line,

def file_len(file):
    with open(file) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def animate(i):
    len = file_len('output/test.txt')
    with open('output/test.txt') as f:
        d=np.genfromtxt(itertools.islice(f,len-50,len),delimiter=',',usecols=(0,1))
    xdata=d[:,0]
    ydata=d[:,1]
    line.set_data(xdata, ydata)
    ax.relim()
    ax.autoscale()
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, interval=25)

plt.show()