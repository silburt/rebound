#This macro calculates the elapsed time versus the number of planetesimals in the simulation. All simulations are integrated for the same amount of time and have the same basic setup. Compare to mercury and swifter.

import glob
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import re
from scipy.optimize import curve_fit

def func(x,a,b):
    return a*x**b

def get_times(files, ext, Np_pos):
    ET = []
    Np = []
    for f in files:
        try:
            lines = open(f+ext, 'r').readlines()
            elapsed = float(lines[0].split()[-2])
            planetesimals = int(f.split("_")[Np_pos].split("Np")[1])
            ET.append(elapsed)
            Np.append(planetesimals)
        except:
            print "couldnt find file",f
    return Np, ET

def get_energy(files, ext):
    dE = []
    if ext == 'hermes':
        for i in range(len(files)):
            files[i] = files[i].split('_ET.txt')[0]+'.txt'
        for f in files:
            lines = open(f, 'r').readlines()
            dE.append(float(lines[-1].split(',')[1]))
    elif ext == 'swifter':
        for f in files:
            lines = open(f+'/energyoutput.txt', 'r').readlines()
            dE.append(float(lines[-1].split()[1]))
    elif ext == 'mercury':
        for f in files:
            lines = open(f+'/eo.txt', 'r').readlines()
            dE.append(float(lines[-2].split()[4]))
    return dE

ms = 0.25
alpha = 0.4

dirP = str(sys.argv[1])

print "...HERMES..."
files = glob.glob(dirP+'*_ET.txt')
Np_H, times_H = get_times(files,'', -3)
dE_H = get_energy(files, 'hermes')

#x,y = zip(*[(x,y) for (x,y) in sorted(zip(Np_H,times_H))])
#popt, pcov = curve_fit(func, x, y)
#plt.plot(x, y, 'o', markeredgecolor='none', color='darkgreen', label='Numerical Simulation')
#plt.plot(x, func(x,popt[0],popt[1]), label='Best Fit: t$_{elapsed}$='+str('%.1e'%popt[0])+'$\cdot N_{pl}$$^{'+str(round(popt[1],2))+'}$')

print "...Swifter..."
dir = '../../../swifter/example/saved_input_files/round7_10yr_allNp_comp/'
files = [x[0] for x in os.walk(dir)][1:]
Np_S, times_S = get_times(files,'/elapsed_time.txt', -2)
dE_S = get_energy(files, 'swifter')

print "...Mercury..."
dir = '../../../mercury6/saved_input_files/round9_10yr_allNp_comp/'
files = [x[0] for x in os.walk(dir)][1:]
Np_M, times_M = get_times(files,'/ET.txt', -2)
dE_M = get_energy(files, 'mercury')

fig = plt.figure(figsize=(8,10))
ax = fig.add_subplot(2, 1, 1)
ax.plot(Np_H, times_H, 'o', markeredgecolor='none', color='darkgreen', label='HERMES')
ax.plot(Np_S, times_S, 'o', markeredgecolor='none', color='darkblue', label='SyMBA')
ax.plot(Np_M, times_M, 'o', markeredgecolor='none', color='darkred', label='MERCURY')
ax.legend(loc='upper left', numpoints=1, fontsize=10)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([0.5,120000])
ax.set_xlabel('Initial Number of Planetesimals',fontsize=13)
ax.set_ylabel('Elapsed Time (seconds)',fontsize=13)

ax = fig.add_subplot(2, 1, 2)
ax.plot(Np_H, dE_H, 'o', markeredgecolor='none', color='darkgreen', label='HERMES')
ax.plot(Np_S, dE_S, 'o', markeredgecolor='none', color='darkblue', label='SyMBA')
ax.plot(Np_M, dE_M, 'o', markeredgecolor='none', color='darkred', label='MERCURY')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim([0.5,120000])
ax.set_xlabel('Initial Number of Planetesimals',fontsize=13)
ax.set_ylabel('Fractional Energy Error',fontsize=13)

plt.savefig(dirP+'ETvNp.pdf')
plt.show()
