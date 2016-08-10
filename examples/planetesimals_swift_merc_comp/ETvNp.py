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

def get_times(files, ext):
    ET = []
    Np = []
    for f in files:
        try:
            ff = open(f+ext, 'r')
            lines = ff.readlines()
            if ext == '':
                elapsed = float(lines[1].split()[-2])/3600.
                planetesimals = int(re.split("_|/",f)[-4].split("Np")[1])
            else:
                elapsed = float(lines[0].split()[-1])/3600.
                #elapsed = (float(lines[1].split()[-1]) - float(lines[0].split()[-1]))/3600.
                planetesimals = int(re.split("_|/",f)[-2].split("Np")[1])
            ET.append(elapsed)
            Np.append(planetesimals)
        except:
            print "couldnt find file",f
    return Np, ET

ms = 0.25
alpha = 0.4

dirP = str(sys.argv[1])

#Elapsed time - HERMES
files = glob.glob(dirP+'*_elapsedtime.txt')
Np_H, times_H = get_times(files,'')
x,y = zip(*[(x,y) for (x,y) in sorted(zip(Np_H,times_H))])
popt, pcov = curve_fit(func, x, y)
plt.plot(x, y, 'o', markeredgecolor='none', color='darkgreen', label='Numerical Simulation')
plt.plot(x, func(x,popt[0],popt[1]), label='Best Fit: t$_{elapsed}$='+str('%.1e'%popt[0])+'$\cdot N_{pl}$$^{'+str(round(popt[1],2))+'}$')

'''
#Elapsed time - Swifter
dir = '../../../swifter/example/input_files/'
files = [x[0] for x in os.walk(dir)][1:]
Np_S, times_S = get_times(files,'/ET.txt')
plt.plot(Np_S, times_S, 'o', markeredgecolor='none', color='darkblue', label='SyMBA')

#Elapsed time - MERCURY
dir = '../../../mercury6/input_files/'
files = [x[0] for x in os.walk(dir)][1:]
Np_M, times_M = get_times(files,'/ET.txt')
plt.plot(Np_M, times_M, 'o', markeredgecolor='none', color='darkred', label='MERCURY')
'''

plt.legend(loc='upper left', numpoints=1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Initial Number of Planetesimals',fontsize=13)
plt.ylabel('Elapsed Time (hours)',fontsize=13)
plt.savefig(dirP+'ETvNp.pdf')
plt.show()