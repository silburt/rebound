#This macro calculates the elapsed time versus the number of planetesimals in the simulation. All simulations are integrated for the same amount of time and have the same basic setup. Compare to mercury and swifter.

import glob
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import re

def get_times(files, ext):
    elapsed = []
    Np = []
    for f in files:
        ff = open(f+ext, 'r')
        lines = ff.readlines()
        if ext == '':
            elapsed.append(float(lines[1].split()[-2])/3600.)
            Np.append(int(re.split("_|/",f)[-4].split("Np")[1]))
        else:
            elapsed.append((float(lines[1].split()[-1]) - float(lines[0].split()[-1]))/3600.)
            Np.append(int(re.split("_|/",f)[-2].split("Np")[1]))
    return Np, elapsed

ms = 0.25
alpha = 0.4

dirP = str(sys.argv[1])

#Elapsed time - HERMES
files = glob.glob(dirP+'*_elapsedtime.txt')
Np_H, times_H = get_times(files,'')
plt.plot(Np_H, times_H, 'o', markeredgecolor='none', color='darkgreen', label='HERMES')


#Elapsed time - Swifter
dir = '../../../swifter/example/input_files/'
files = [x[0] for x in os.walk(dir)][1:]
Np_S, times_S = get_times(files,'/elapsed_time.txt')
plt.plot(Np_S, times_S, 'o', markeredgecolor='none', color='darkblue', label='SyMBA')

#Elapsed time - MERCURY
dir = '../../../mercury6/input_files/'
files = [x[0] for x in os.walk(dir)][1:]
Np_M, times_M = get_times(files,'/elapsed_time.txt')
plt.plot(Np_M, times_M, 'o', markeredgecolor='none', color='darkred', label='MERCURY')

plt.xscale('log')
plt.xlabel('Number of Planetesimals',fontsize=13)
plt.ylabel('Elapsed Time (hours)',fontsize=13)
plt.show()