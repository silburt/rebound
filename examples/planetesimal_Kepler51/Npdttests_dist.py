#The purpose of this macro is to show the distribution of systems in Np/dt space, and color code according to their period ratios.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob
import sys

dir = sys.argv[1]
files = glob.glob(dir+'*info.txt')
N = len(files)

Np = np.zeros(N)
dt = np.zeros(N)
Pratio = np.zeros(N)

for i in range(N):
    try:
        f = files[i]
        Np[i] = float(f.split("_")[-3].split("Np")[1])
        dt[i] = float(f.split("_")[-4].split("dt")[1])
        fos = open(files[i].split("_info.txt")[0]+".txt", 'r')
        time, dE, N, N_mini, HSF, m1, a1, e1, m2, a2, e2, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)
        Pratio[i] = (np.mean(a2[-1]/a1[-1]))**(1.5)
    except:
        Np[:-1]
        dt[:-1]
        Pratio[:-1]
        print files[i], "had an ejected particle or something"

plt.scatter(Np, dt, c=Pratio, cmap = cm.rainbow, lw=0)
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-5,5e-2])
plt.xlim([1e2,1e5])
plt.xlabel('Number of planetesimals')
plt.ylabel('timestep (years)')
plt.colorbar()
plt.show()
