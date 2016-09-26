from matplotlib import *
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

res_thresh = 0.05

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

data  = np.genfromtxt("planets.csv", dtype=None, delimiter=',',comments='#',filling_values=0)
#data = data[data[:,3]=="Transit"]

row, name, letter, method, Npl, P, Pup, Pdown, junk1, a, aup, adown, junk2, e, eup, edown, junk3, m, mup, mdown, junk4, mprov, r, rup, rdown, junk5, ttvflag, d, dup, ddown, junk6, Ms, Msup, Msdown, junk7, junk8, Rs, Rsup, Rsdown, junk9, junk10 = zip(*data)

name_current = ''
i=0
while i < len(data):
    if name_current != name[i]:
        name_current = name[i]
        Parr, marr, rarr = np.asarray(zip(*sorted(zip(P[i:i+Npl[i]],m[i:i+Npl[i]],r[i:i+Npl[i]]))))
        Pratio = Parr[1:-1]/Parr[0:-2]
        index = np.where(np.abs(Pratio - 2.0) < res_thresh)[0]
        index = np.concatenate((index, np.where(np.abs(Pratio - 1.5) < res_thresh)[0]))
        for idx in index:
            ax.scatter(Pratio[idx],rarr[idx],rarr[idx+1])
            ax.text(Pratio[idx],rarr[idx],rarr[idx+1],'%s'%(name[i]), size=10, zorder=1, color='k')
        i += Npl[i]
    else:
        i += 1
        print "Something weird happened here"

ax.set_xlabel('Period Ratio')
ax.set_ylabel('Radius Planet 1 (R$_{Jupiter}$)')
ax.set_zlabel('Radius Planet 2 (R$_{Jupiter}$)')
ax.set_ylim([0,0.5])
ax.set_zlim([0,0.5])

plt.show()