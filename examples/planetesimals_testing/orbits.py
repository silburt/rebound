import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

diagnostics = 1
plot_dr_or_com = 1  #0 = dr, 1 = com

fos = open(''+file_name, 'r')
data = np.loadtxt(fos, delimiter=',')

ms=3
if diagnostics == 1:
    time = data[:,0]
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(10,10))
    axes[0].plot(time,data[:,1], 'o', ms=ms, markeredgecolor='none')
    axes[0].plot(time,0.5e-12*time, color='red', label='t')
    axes[0].plot(time,1e-12*time**0.5, color='black', label='t^0.5')
    axes[0].legend(loc='upper left',prop={'size':10})
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    axes[0].set_xlim([0.1,max(time)])
    axes[0].set_ylabel('Energy')
    axes[1].plot(time,data[:,2], 'o', ms=ms, markeredgecolor='none', label='global: Np')
    axes[1].plot(time,data[:,3], 'o', ms=ms, markeredgecolor='none', label='mini: Np')
    axes[1].legend(loc='upper left',prop={'size':10})
    axes[1].set_xscale('log')
    axes[1].set_ylim([0,max(data[:,2])])
    axes[1].set_ylabel('Number of particles')
    #axes[2].plot(time,data[:,10], 'o', ms=ms, markeredgecolor='none')
    #axes[2].set_ylabel('arg of peri single particle')
    #axes[2].set_ylim([0,6.28])
    #axes[2].plot(time,data[:,4], 'o', ms=ms, markeredgecolor='none',label='elapsed time (s)')
    #axes[2].plot(time,data[:,5], 'o', ms=ms, markeredgecolor='none',label='N Close Encounters')
    #axes[2].legend(loc='upper left',prop={'size':10})
    #axes[2].set_xscale('log')
    axes[2].plot(time,data[:,11], 'o', ms=ms, markeredgecolor='none',label='')
    axes[2].plot(time,1e-11*time**1.5, color='black',label='t^1.5')
    axes[2].plot(time,1e-11*time**2, color='red',label='t^2')
    axes[2].set_xlabel('simulation time (yrs)')
    axes[2].set_ylabel('dcom')
    axes[2].set_yscale('log')
    axes[2].legend(loc='lower left',prop={'size':10})
else:
    plt.plot(data[:,0],data[:,1], 'o', ms=ms, markeredgecolor='none')
    #plt.plot(data[:,0],3e-10*data[:,0]**(0.5),color='black',label='t^1/2 growth')
    plt.ylabel('Energy')
    plt.xlabel('time (years)')
    #plt.legend(loc='upper left',prop={'size':10})
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([0.1,data[-1,0]])

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
