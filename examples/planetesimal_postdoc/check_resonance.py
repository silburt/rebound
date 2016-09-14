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
time, dE, N, N_mini, a1, e1, m1, a2, e2, m2, a3, e3, m3, a4, e4, m4, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)

ms=3
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(10,10))
plt.subplots_adjust(hspace = 0.35)

#periods
P1 = np.sqrt(4*np.pi**2 * a3**3)
P2 = np.sqrt(4*np.pi**2 * a4**3)

#plot
axes[0].plot(time,P2/P1,'o',ms=ms, markeredgecolor='none')
axes[1].plot(time, e1, 'o', ms=ms, markeredgecolor='none', label='Earth')
axes[1].plot(time, e2, 'o', ms=ms, markeredgecolor='none', label='Mars')
axes[1].plot(time, e3, 'o', ms=ms, markeredgecolor='none', label='Jupiter')
axes[1].plot(time, e4, 'o', ms=ms, markeredgecolor='none', label='Saturn')
axes[2].plot(time,phi1, 'o', ms=ms, markeredgecolor='none')
axes[3].plot(time,phi2, 'o', ms=ms, markeredgecolor='none')

#labelling
axes[0].set_ylabel('Period Ratio', fontsize=13)
axes[0].set_xlabel('Time (Years)', fontsize=13)
axes[1].set_ylabel('Eccentricity', fontsize=13)
axes[1].set_xlabel('Time (Years)', fontsize=13)
axes[1].legend(loc='upper right')
axes[2].set_ylabel('Phi1', fontsize=13)
axes[2].set_xlabel('Time (Years)', fontsize=13)
axes[3].set_ylabel('Phi2', fontsize=13)
axes[3].set_xlabel('Time (Years)', fontsize=13)

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'_rescheck.png')
plt.show()
