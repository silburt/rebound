from matplotlib import *
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111, projection='3d')

#data  = np.genfromtxt("planets.csv", dtype=None, delimiter=',',comments='#',filling_values=0)
#row, name, letter, method, Npl, P, Pup, Pdown, junk1, a, aup, adown, junk2, e, eup, edown, junk3, m, mup, mdown, junk4, mprov, r, rup, rdown, junk5, ttvflag, d, dup, ddown, junk6, Ms, Msup, Msdown, junk7, junk8, Rs, Rsup, Rsdown, junk9, junk10 = np.asarray(zip(*data))
header = ['name', 'letter', 'method', 'Npl', 'P', 'Pup', 'Pdown', 'a', 'aup', 'adown', 'e', 'eup', 'edown', 'm', 'mup', 'mdown', 'r', 'rup', 'rdown', 'Ms', 'Msup', 'Msdown', 'Rs', 'Rsup', 'Rsdown']
usecols = [1,2,3,4,5,6,7,9,10,11,13,14,15,17,18,19,22,23,24,31,32,33,36,37,38]
data = pd.read_csv("planets.csv",names=header, delimiter=',',usecols=usecols,skiprows=48)

#FOR ML restricting to 3-planet systems
#data = data[data["Npl"] == 3]

#Parameters
var = 'm'           #variable to plot
res_thresh = 1      #distance from resonance.
res = 2

systems = np.unique(data['name'])

for n in systems:
    index = data['name']==n
    Parr = data.loc[index,'P'].values
    vararr = data.loc[index,var].values
    Pratio = Parr[1:]/Parr[0:-1]
    index2 = np.where(np.abs(Pratio - res) < res_thresh)[0]
    for idx in index2:
        #if vararr[idx] > 0 and vararr[idx+1] > 0:
        if vararr[idx] > 0.5 and vararr[idx+1] > 0.5:
            if var == 'm':
                ax.scatter(Pratio[idx],np.log10(vararr[idx]), np.log10(vararr[idx+1]))
                ax.text(Pratio[idx],np.log10(vararr[idx]), np.log10(vararr[idx+1]),'%s'%(n), size=10, zorder=1, color='k')
            else:
                ax.scatter(Pratio[idx],vararr[idx], vararr[idx+1])
                ax.text(Pratio[idx],vararr[idx], vararr[idx+1],'%s'%(n), size=7, zorder=1, color='k')

if var == 'm':
    ax.set_ylabel('Mass Planet 1 (log10(M$_{Jupiter}$))')
    ax.set_zlabel('Mass Planet 2 (log10(M$_{Jupiter}$))')
elif var == 'r':
    ax.set_ylabel('Radius Planet 1 (R$_{Jupiter}$)')
    ax.set_zlabel('Radius Planet 2 (R$_{Jupiter}$)')

ax.set_xlabel('Period Ratio')
ax.set_title('Planets with period ratio close to: %.3f'%res)
plt.show()
