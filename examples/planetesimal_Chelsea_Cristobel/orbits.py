import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

file_name=str(sys.argv[1])
fos = open(file_name, 'r')
data = np.loadtxt(fos, delimiter=',')

plt.plot(data[:,0],data[:,1], 'o', markeredgecolor='none')
plt.ylabel('Energy')
plt.xlabel('time (years)')
#plt.legend(loc='upper left',prop={'size':10})
plt.xscale('log')
plt.yscale('log')
plt.xlim([0.1,data[-1,0]])
plt.show()
