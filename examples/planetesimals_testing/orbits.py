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

fos = open(file_name, 'r')
data = np.loadtxt(fos, delimiter=',')

signed_energy = 1
reverse_xaxis = 0

ms=3
if diagnostics == 1:
    time,dE = data[:,0],data[:,1]
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(10,10))
    if signed_energy == 1:
        axes[0].plot(time[dE>0],abs(dE[dE>0]),'o', ms=ms, markeredgecolor='none',color='blue',label='positive energy')
        axes[0].plot(time[dE<0],abs(dE[dE<0]),'o', ms=ms, markeredgecolor='none',color='green',label='negative energy')
        energytitle = 'signed fractional energy'
    else:
        axes[0].plot(time,abs(dE), 'o', ms=ms, markeredgecolor='none',color='blue')
        axes[0].plot(time,0.5e-12*time, color='red', label='t')
        axes[0].plot(time,1e-12*time**0.5, color='black', label='t^0.5')
        energytitle = 'fractional energy'
        axes[0].set_xscale('log')
    axes[0].set_xlim([0.1,max(time)])
    if reverse_xaxis == 1:
        axes[0].set_xlim([max(time),0.1])
    axes[0].set_yscale('log')
    axes[0].legend(loc='upper left',prop={'size':10})
    axes[0].set_ylabel(energytitle)
    axes[1].plot(time,data[:,2], 'o', ms=ms, markeredgecolor='none', label='global: Np')
    axes[1].plot(time,data[:,3], 'o', ms=ms, markeredgecolor='none', label='mini: Np')
    axes[1].legend(loc='upper left',prop={'size':10})
    axes[1].set_yscale('log')
    axes[1].set_ylim([0.1,max(data[:,2])])
    axes[1].set_ylabel('Number of particles')
#HSR
    axes[2].plot(time,data[:,10], 'o', ms=ms, markeredgecolor='none')
    axes[2].set_ylabel('HSR')
#com
    #axes[2].plot(time,data[:,11], 'o', ms=ms, markeredgecolor='none',label='')
    #axes[2].plot(time,1e-11*time**1.5, color='black',label='t^1.5')
    #axes[2].plot(time,1e-11*time**2, color='red',label='t^2')
    #axes[2].set_xlabel('simulation time (yrs)')
    #axes[2].set_ylabel('dcom')
    #axes[2].set_yscale('log')
    #axes[2].legend(loc='lower left',prop={'size':10})
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
