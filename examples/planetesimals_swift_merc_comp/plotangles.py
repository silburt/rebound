#plots the distribution of incoming/outgoing angles to the planet's hill sphere.

import matplotlib.pyplot as plt
import numpy as np
import sys

def sort_index(data):
    N = len(data)
    sorted = np.zeros(N);
    for i in xrange(0,N):
        index = np.argmin(data) #index of min
        sorted[i] = data[index]*57.29578
        data = np.delete(data,index)
    y = np.arange(1,N+1)/float(N)
    return sorted, y

filename=str(sys.argv[1])
filename = filename.split(".txt")
fos = open(filename[0]+'_inangle.txt', 'r')
datain = np.loadtxt(fos, delimiter=',')
fos = open(filename[0]+'_outangle.txt', 'r')
dataout = np.loadtxt(fos, delimiter=',')

xin, yin = sort_index(datain[:,1])
xout, yout = sort_index(dataout[:,1])

fig, axes = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(8,8))
axes[0].plot(xin,yin)
axes[1].plot(xout,yout)
axes[1].set_xlabel('outgoing planetesimal angle (degrees)')
axes[0].set_xlabel('incoming planetesimal angle (degrees)')
axes[0].set_xlim([0,360])
axes[0].set_ylim([0,1])
axes[1].set_ylim([0,1])
print 'Preparing PDF'
plt.savefig(filename[0]+'_plotangles.png')
plt.show()