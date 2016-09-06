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

Npl_cutoff = 100000

ms = 0.25
alpha = 0.4

dir = sys.argv[1]
files = glob.glob(dir+'*sd*.txt')
N = len(files)
i=0
while i < N:    #just want the main .txt files
    f = files[i]
    string = f.split("_")
    try:
        Npl = float(string[-2].split("Np")[1])
    except:
        Npl = 10e10
    if string[-1]=="info.txt" or string[-1]=="elapsedtime.txt" or string[-2]=="eiasnapshot" or Npl>Npl_cutoff:
        files.remove(files[i])
        N -= 1
    else:
        i += 1

files = sorted(files, key = natural_key)
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10,10), sharex=True)

for f in files[5:]:
    time, dE, N, N_mini, a1, e1, a2, e2, phi1, phi2, phi3 = np.loadtxt(f, delimiter=',', unpack=True)
    axes[0].plot(time,np.sqrt((a2**3)/(a1**3)), label='N='+str(N[0]))
    axes[1].plot(time,N_mini/N[0],'.')
axes[0].set_ylabel('P2/P1')
axes[0].set_xlim([3e3,1e6])
axes[0].legend(loc='upper left',numpoints=1)
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[1].set_xlabel('time (years)')
axes[1].set_ylabel('$N_{mini}/N$')
axes[1].set_ylim([1e-4,0.5])
plt.savefig(dir+'PvNp.png')

#plt.show()