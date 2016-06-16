import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

fos = open(''+file_name, 'r')
time, dE, N, N_mini, a1, e1, a2, e2, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)
time /= 6.2831853

ms=3
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(10,10))
plt.subplots_adjust(hspace = 0.25)

#periods
P1 = np.sqrt(4*np.pi**2 * a1**3)
P2 = np.sqrt(4*np.pi**2 * a2**3)

#plot
axes[0].plot(time, dE, 'o', color='blue', ms=ms, markeredgecolor='none')
axes[1].plot(time,P2/P1,'o',ms=ms, markeredgecolor='none')
axes[2].scatter(e1*np.cos(phi1),e1*np.sin(phi1),c=time, cmap=cm.rainbow)
axes[3].scatter(e1*np.cos(phi2),e1*np.sin(phi2),c=time, cmap=cm.rainbow)

#labelling
axes[0].set_ylabel('Fractional Energy Error', fontsize=13)
axes[0].set_yscale('log')
axes[0].set_ylabel('Fractional', fontsize=13)
axes[0].set_xlabel('Time (Years)', fontsize=13)
axes[1].set_ylabel('Period Ratio (AU)', fontsize=13)
axes[1].set_xlabel('Time (Years)', fontsize=13)
axes[2].set_ylabel('e1*cos(phi1)', fontsize=13)
axes[2].set_xlabel('e1*sin(phi1)', fontsize=13)
axes[3].set_ylabel('e1*cos(phi2)', fontsize=13)
axes[3].set_xlabel('e1*sin(phi2)', fontsize=13)

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
