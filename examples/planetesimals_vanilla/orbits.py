import sys
import matplotlib.pyplot as plt
import numpy as np
import re

file_name = "output/test.txt"
fos = open(file_name, 'r')
time, dE, ET, HSF, N, Nmini, miniactivesteps = np.loadtxt(fos, delimiter=',', unpack=True)
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,10), sharex=True)

axes[0].plot(time,dE, '.')
axes[0].plot(time,0.5e-11*time, color='red', label='t')
axes[0].plot(time,1e-10*time**0.5, color='black', label='t^0.5')
axes[0].set_xscale('log')
axes[0].set_xlim([0.1,max(time)])
axes[0].set_yscale('log')
axes[1].plot(time,miniactivesteps, '.')
axes[1].set_ylabel('Total steps mini active')
axes[1].set_yscale('log')
axes[2].plot(time,ET)
axes[2].set_ylabel('Elapsed Time')
axes[2].set_yscale('log')

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
