#This macro calculates the average energy for a set of HERMES runs.
#This plot is going in the paper

import glob
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import re
from scipy.optimize import curve_fit

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def line(x,m,b):
    return m*x + b

ms = 0.25
alpha = 0.4
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10,10))

#Elapsed time
dirP = str(sys.argv[1])
files = glob.glob(dirP+'*_elapsedtime.txt')
files = sorted(files, key = natural_key)
elapsed = []
Np = []
for f in files:
    ff = open(f, 'r')
    lines = ff.readlines()
    elapsed.append(float(lines[1].split()[-2])/3600.)
    Np.append(float(re.split("_|/",f)[-4].split("Np")[1]))
Np = np.asarray(Np, dtype='float')
elapsed = np.asarray(elapsed)

logNp = np.log10(Np)
logelapsed = np.log10(elapsed)
popt, pcov = curve_fit(line, logNp, logelapsed) #find the fit in log space

axes[0].plot(Np,elapsed, 'o', label='Simulation')
axes[0].plot(Np,10**(popt[1])*Np**popt[0],label='t='+str("%.2e"%10**popt[1])+'$N_p^{'+str(round(popt[0],2))+'}$')
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[0].set_xlabel('Number of Planetesimals',fontsize=13)
axes[0].set_ylabel('Elapsed Time (hours)',fontsize=13)
axes[0].legend(loc='upper left',numpoints=1)
axes[0].plot(Np,0.4e-3*Np)

N_plot = np.array([13,30,82,145,255,449,1048,3237,6551,8685]) + 2
for f in files:
    time, dE, a, N, N_mini, mini_on, ET = np.loadtxt(f.split("_elapsedtime.txt")[0]+".txt", delimiter=',', unpack=True)
    if N[0] in N_plot:
        axes[1].plot(time,N_mini/N[0],'.', label='N='+str(N[0]))
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[1].set_xlabel('time (years)')
axes[1].set_ylabel('$N_{mini}/N$')
axes[1].legend(loc='upper left',numpoints=1)

plt.savefig(dirP+'ETvNp.png')
plt.show()