from matplotlib import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

data  = np.genfromtxt("planets.csv", dtype=None, delimiter=',',comments='#',filling_values=0)
row, name, letter, method, Npl, P, Pup, Pdown, junk1, a, aup, adown, junk2, e, eup, edown, junk3, m, mup, mdown, junk4, mprov, r, rup, rdown, junk5, ttvflag, d, dup, ddown, junk6, Ms, Msup, Msdown, junk7, junk8, Rs, Rsup, Rsdown, junk9, junk10 = np.asarray(zip(*data))

m=m.astype(np.float)
r=r.astype(np.float)
P=P.astype(np.float)

#Parameters
var = m             #variable to plot
res_thresh = 0.06   #distance from resonance.
res = 2

systems = np.unique(name)

for n in systems:
    index = np.where(name==n)
    Parr = P[index]
    vararr = var[index]
    Pratio = Parr[1:]/Parr[0:-1]
    index2 = np.where(np.abs(Pratio - res) < res_thresh)[0]
    for idx in index2:
        if vararr[idx] > 0 and vararr[idx+1] > 0:
            ax.scatter(Pratio[idx],np.log10(vararr[idx]), np.log10(vararr[idx+1]))
            ax.text(Pratio[idx],np.log10(vararr[idx]), np.log10(vararr[idx+1]),'%s'%(n), size=10, zorder=1, color='k')

if var.all() == m.all():
    ax.set_ylabel('Mass Planet 1 (log10(M$_{Jupiter}$))')
    ax.set_zlabel('Mass Planet 2 (log10(M$_{Jupiter}$))')
elif var.all() == r.all():
    ax.set_ylabel('Radius Planet 1 (R$_{Jupiter}$)')
    ax.set_zlabel('Radius Planet 2 (R$_{Jupiter}$)')
    ax.set_ylim([0,0.5])
    ax.set_zlim([0,0.5])

ax.set_xlabel('Period Ratio')
ax.set_title('Planets with period ratio close to: %.3f'%res)
plt.show()