import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

file_name=str(sys.argv[1])
fos = open(file_name, 'r')
data = np.loadtxt(fos, delimiter=',')

fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(10,10))

axes[0].plot(data[:,0],data[:,1], 'o', markeredgecolor='none')
axes[0].set_ylabel('Energy')
axes[0].set_xlabel('time (years)')
#axes[0].legend(loc='upper left',prop={'size':10})
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[0].set_xlim([0.1,data[-1,0]])
axes[1].plot(data[:,0],data[:,4],'ko',markeredgecolor='none',label='N_tot')
axes[1].plot(data[:,0],data[:,2],'o',markeredgecolor='none',label='N_mini')
axes[1].legend(loc='upper left',prop={'size':10})
axes[1].set_ylim([0,10])
axes[2].plot(data[:,0],data[:,6],label='time')
axes[2].set_ylabel('elapsed time (s)')

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
