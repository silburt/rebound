import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

diagnostics = 1

fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=',')

ms=3
if diagnostics == 1:
    time = data[:,0]
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,10))
    axes[0].plot(time,data[:,1], 'o', ms=ms, markeredgecolor='none')
    axes[0].plot(time,1e-7*time, color='red', label='t')
    axes[0].plot(time,1e-6*time**0.5, color='black', label='t^0.5')
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    axes[0].set_xlim([1,max(time)])
    axes[0].set_ylabel('dE/E$_0$')
    axes[1].plot(time,data[:,2], 'o', ms=ms, markeredgecolor='none', label='global: Np')
    axes[1].plot(time,data[:,3]-3, 'o', ms=ms, markeredgecolor='none', label='mini: Np')
    axes[1].legend(loc='upper left',prop={'size':10})
    axes[1].set_xscale('log')
    axes[1].set_ylim([0,max(data[:,2])])
    axes[1].set_xlim([1,max(time)])
    axes[1].set_ylabel('Number of particles')
    axes[1].set_xscale('log')
    axes[2].plot(time,data[:,7], 'o', ms=ms, markeredgecolor='none')
    axes[2].set_ylabel('Nsteps Mini Active')
    axes[2].set_yscale('log')
    axes[2].set_xscale('log')
    print "percent mini active is %f"%(sum(data[:,4])/len(data[:,4]))
else:
    plt.plot(data[:,0],data[:,1], 'o', ms=ms, markeredgecolor='none')
    plt.plot(data[:,0],3e-10*data[:,0]**(0.5),color='black',label='t^1/2 growth')
    plt.ylabel('dE/E$_0$')
    plt.xlabel('time (years)')
    plt.legend(loc='upper left',prop={'size':10})
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([0.5,data[-1,0]])

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
#plt.show()
