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
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(10,10), sharex=True)
plt.subplots_adjust(hspace = 0.35)

#plot
axes[0].plot(time, dE,'o',ms=ms, markeredgecolor='none')
axes[0].plot(time, 1e-8*time**0.5)
axes[1].plot(time, N, 'o', ms=ms, markeredgecolor='none')
axes[1].plot(time, N_mini, 'o', ms=ms, markeredgecolor='none')
axes[2].plot(time, HSF, 'o', ms=ms, markeredgecolor='none')
axes[3].plot(time, (a1 - a1.mean())/a1.std(), 'o', ms=ms, markeredgecolor='none',label='planet 1')
axes[3].plot(time, (a2 - a2.mean())/a2.std(), 'o', ms=ms, markeredgecolor='none',label='planet 2')
if Nplanets == 3:
    axes[3].plot(time, (a3 - a3.mean())/a3.std(), 'o', ms=ms, markeredgecolor='none',label='planet 3')
#axes[3].plot(time, e1, 'o', ms=ms, markeredgecolor='none', label='planet 1')
#axes[3].plot(time, e2, 'o', ms=ms, markeredgecolor='none', label='planet 2')
#axes[3].plot(time, e3, 'o', ms=ms, markeredgecolor='none', label='planet 3')

#labelling
axes[0].set_ylabel('Fractional Energy', fontsize=13)
axes[0].set_yscale('log')
axes[0].set_xscale('log')
axes[0].set_xlim([100,time[-1]])
axes[1].set_ylabel('Number of particles', fontsize=13)
axes[1].set_yscale('log')
axes[1].set_xlabel('Time (Years)', fontsize=13)
axes[2].set_ylabel('Hybrid Switch Ratio', fontsize=13)
axes[2].set_xlabel('Time (Years)', fontsize=13)
axes[3].set_ylabel(r'$(a-\bar{a})/\sigma_a$', fontsize=13)
axes[3].legend(loc='upper left')
axes[3].set_xlabel('Time (Years)', fontsize=13)
#axes[3].set_ylabel('Eccentricity', fontsize=13)
#axes[3].set_xlabel('Time (Years)', fontsize=13)

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'_orbit.png')
#plt.show()
