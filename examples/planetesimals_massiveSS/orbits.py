import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

fos = open(''+file_name, 'r')
time, dE, N, a1, a2, a3, a4 = np.loadtxt(fos, delimiter=',', unpack=True)
color = ['yellow','red','green','blue','purple']
labels = ['Jupiter','Saturn','Uranus','Neptune']

ms=3
fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(8,6))
plt.subplots_adjust(hspace = 0.1)
t=0
i=0

time /= 6.2831853
axes[0].plot(time, dE, 'o', color='blue', ms=ms, markeredgecolor='none')

#stitch together arrays that got screwed up.

i1 = np.where(a4!=0)
i2 = np.where(a4==0)
jupiter = a1
saturn = a2[i1]
uranus = np.concatenate((a3[i1],a2[i2]))
neptune = np.concatenate((a4[i1],a3[i2]))

#plot
axes[1].plot(time,jupiter,'o',ms=ms, markeredgecolor='none',color=color[0],label=labels[0])
axes[1].plot(time[i1],saturn,'o',ms=ms, markeredgecolor='none',color=color[1],label=labels[1])
axes[1].plot(time,uranus,'o',ms=ms, markeredgecolor='none',color=color[2],label=labels[2])
axes[1].plot(time,neptune,'o',ms=ms, markeredgecolor='none',color=color[3],label=labels[3])

axes[0].set_ylabel('Fractional Energy Error', fontsize=13)
axes[0].set_yscale('log')
axes[1].legend(loc='upper right',prop={'size':10}, numpoints=1, markerscale=3)
axes[1].set_ylabel('Semi-Major Axis (AU)', fontsize=13)
axes[1].set_xlabel('Time (Years)', fontsize=13)
axes[1].set_ylim([0,50])

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
