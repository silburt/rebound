import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

fos = open(''+file_name, 'r')
time, dE, N, junk, junk2, a1, a2, a3, a4 = np.loadtxt(fos, delimiter=',', unpack=True)
color = ['yellow','red','green','blue','purple']
labels = ['Jupiter','Saturn','Uranus','Neptune']

ms=3
fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(8,8))
plt.subplots_adjust(hspace = 0.1)
t=0
i=0

a = [a1,a2,a3,a4]
time /= 6.2831853
axes[0].plot(time, dE, 'o', color='blue', ms=ms, markeredgecolor='none')
for i in xrange(0,len(labels)):
    axes[1].plot(time,a[i],'o',ms=ms, markeredgecolor='none',color=color[i],label=labels[i])

axes[0].set_ylabel('dE/E(0)', fontsize=16)
axes[0].set_yscale('log')
axes[1].legend(loc='upper right',prop={'size':10}, numpoints=1, markerscale=3)
axes[1].set_ylabel('semi-major axis', fontsize=16)
axes[1].set_xlabel('time (years)', fontsize=16)
axes[1].set_xlim([0,1000])
axes[1].set_ylim([0,120])

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
