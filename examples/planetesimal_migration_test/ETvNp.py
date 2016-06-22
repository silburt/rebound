#This macro calculates the average energy for a set of PLANETESIMAL, SWIFTER and MERCURY runs and plots them all on one sexy figure.
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

plt.plot(Np, elapsed, 'o')
plt.xscale('log')
plt.xlabel('Number of Planetesimals',fontsize=13)
plt.ylabel('Elapsed Time (hours)',fontsize=13)
plt.show()