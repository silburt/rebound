#This macro calculates the average energy for a set of HERMES runs.
#This plot is going in the paper

import glob
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import re

ms = 0.25
alpha = 0.4

#Elapsed time
dirP = str(sys.argv[1])
files = glob.glob(dirP+'*_elapsedtime.txt')
elapsed = []
Np = []
for f in files:
    ff = open(f, 'r')
    lines = ff.readlines()
    elapsed.append(float(lines[1].split()[-2])/3600.)
    Np.append(float(re.split("_|/",f)[-4].split("Np")[1]))

plt.plot(Np, elapsed, 'o', label='HSF=3, dt=0.05')
plt.plot(Np, 3e-4*np.asarray(Np)**0.8, label='3e-4*Np^0.8')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Number of Planetesimals',fontsize=13)
plt.ylabel('Elapsed Time (hours)',fontsize=13)
plt.legend(loc='upper left',numpoints=1)
plt.savefig(dirP+'ETvNp.png')
plt.show()