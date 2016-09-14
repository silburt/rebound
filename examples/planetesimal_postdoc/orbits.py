#This checks whether the planets are still in resonance or not. 
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

fos = open(''+file_name, 'r')
N_bodies = 7
names = ['Earth','Mars','super-Earth','Jupiter','Saturn','Uranus','Neptune']
data = np.loadtxt(fos, delimiter=',')
time, dE, N, N_mini = data[:,0], data[:,1], data[:,2], data[:,3]

ms=3
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10,10), sharex=True)
plt.subplots_adjust(hspace = 0.35)

for i in range(N_bodies):
    axes[0].plot(time,np.sqrt(data[:,3*i+4]**3),'o',ms=ms, markeredgecolor='none',label=names[i])
    axes[1].plot(time,data[:,3*i+5], 'o', ms=ms, markeredgecolor='none', label=names[i])
    #axes[2].plot(time,data[:,3*i+6], 'o', ms=ms, markeredgecolor='none', label=names[i])

axes[0].set_xlim([0,time[-2]])
axes[0].plot([0,max(time)], [12,12], ms=ms, color='black', label='Current Orbit of Jupiter')
axes[0].plot([0,max(time)], [30,30], 'k--', label='Current Orbit of Saturn')
axes[0].set_ylabel('Period (years)', fontsize=13)
axes[0].legend(loc='upper left', fontsize=8, numpoints=1)
axes[1].set_ylabel('eccentricity', fontsize=13)
axes[1].legend(loc='upper left', fontsize=8, numpoints=1)
axes[1].set_yscale('log')
axes[1].set_ylim([1e-3,1])
#axes[2].set_ylabel('mass', fontsize=13)
#axes[2].legend(loc='upper left', fontsize=8, numpoints=1)
#axes[2].set_xlabel('Time (Years)', fontsize=13)


file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'_orbit.png')
#plt.show()