#This checks whether the planets are still in resonance or not. 
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])
Nplanets=int(sys.argv[2])
skip = 0

fos = open(file_name, 'r')
if Nplanets == 2:
    try:
        time, dE, N, N_mini, HSF, m1, a1, e1, m2, a2, e2, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)
    except:
        skip = 1
        print "ejected planet for system %s"%file_name
else:
    time, dE, N, N_mini, HSF, m1, a1, e1, m2, a2, e2, m3, a3, e3, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)

ms=3
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(10,10))
plt.subplots_adjust(hspace = 0.35)

#periods - resonance is between middle/outer planet
if Nplanets == 2:
    Pratio = (a2/a1)**(1.5)
else:
    Pratio = (a3/a2)**(1.5)

#plot
axes[0].plot(time,Pratio,'o',ms=ms, markeredgecolor='none')
axes[1].plot(time, e1, 'o', ms=ms, markeredgecolor='none', label='inner')
axes[1].plot(time, e2, 'o', ms=ms, markeredgecolor='none', label='middle')
if Nplanets == 3:
    axes[1].plot(time, e3, 'o', ms=ms, markeredgecolor='none', label='outer')
axes[2].plot(time,phi1, 'o', ms=ms, markeredgecolor='none')
axes[3].plot(time,phi2, 'o', ms=ms, markeredgecolor='none')
#im = axes[2].scatter(e1*np.cos(phi1),e1*np.sin(phi1),c=time, cmap=cm.rainbow,lw=0)
#axes[3].scatter(e1*np.cos(phi2),e1*np.sin(phi2),c=time, cmap=cm.rainbow,lw=0)
#plt.colorbar(im, label='elapsed time (yr)')

#labelling
axes[0].set_ylabel('Period Ratio', fontsize=13)
axes[0].set_xlabel('Time (Years)', fontsize=13)
axes[1].set_ylabel('Eccentricity', fontsize=13)
axes[1].set_xlabel('Time (Years)', fontsize=13)
axes[1].legend(loc='upper right', numpoints=1, fontsize=8)
axes[2].set_ylabel('Phi1', fontsize=13)
axes[2].set_xlabel('Time (Years)', fontsize=13)
axes[3].set_ylabel('Phi2', fontsize=13)
axes[3].set_xlabel('Time (Years)', fontsize=13)
#axes[2].set_ylabel('e1*cos(phi1)', fontsize=13)
#axes[2].set_xlabel('e1*sin(phi1)', fontsize=13)
#axes[3].set_ylabel('e1*cos(phi2)', fontsize=13)
#axes[3].set_xlabel('e1*sin(phi2)', fontsize=13)

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'_rescheck.png')
#plt.show()
