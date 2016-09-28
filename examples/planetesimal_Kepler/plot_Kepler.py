from matplotlib import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

data  = np.genfromtxt("planets.csv", dtype=None, delimiter=',',comments='#',filling_values=0)
row, name, letter, method, Npl, P, Pup, Pdown, junk1, a, aup, adown, junk2, e, eup, edown, junk3, m, mup, mdown, junk4, mprov, r, rup, rdown, junk5, ttvflag, d, dup, ddown, junk6, Ms, Msup, Msdown, junk7, junk8, Rs, Rsup, Rsdown, junk9, junk10 = zip(*data)

#Parameters
var = m             #variable to plot
res_thresh = 0.05   #distance from resonance.

#analysis
name_current = ''
i=0
while i < len(data):
    if name_current != name[i]:
        name_current = name[i]
        N = 0
        try:
            while name[i+N] == name_current:
                N += 1
        except:
            N = Npl[i]
        Parr, vararr = np.asarray(zip(*sorted(zip(P[i:i+N],var[i:i+N]))))
        Pratio = Parr[1:-1]/Parr[0:-2]
        index = np.where(np.abs(Pratio - 2.0) < res_thresh)[0]
        index = np.concatenate((index, np.where(np.abs(Pratio - 1.5) < res_thresh)[0]))
        for idx in index:
            if vararr[idx] > 0 and vararr[idx+1] > 0:
                ax.scatter(Pratio[idx],np.log10(vararr[idx]), np.log10(vararr[idx+1]))
                ax.text(Pratio[idx],np.log10(vararr[idx]), np.log10(vararr[idx+1]),'%s'%(name[i]), size=10, zorder=1, color='k')
        i += N
    else:
        i += 1
        print "Something weird happened here"

if var == m:
    ax.set_ylabel('Mass Planet 1 (log10(M$_{Jupiter}$))')
    ax.set_zlabel('Mass Planet 2 (log10(M$_{Jupiter}$))')
elif var == r:
    ax.set_ylabel('Radius Planet 1 (R$_{Jupiter}$)')
    ax.set_zlabel('Radius Planet 2 (R$_{Jupiter}$)')
    ax.set_ylim([0,0.5])
    ax.set_zlim([0,0.5])

ax.set_xlabel('Period Ratio')
plt.show()