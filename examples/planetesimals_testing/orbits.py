import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

def cdf(array):
    return np.arange(1,len(array)+1)/float(len(array))

def sort_index(array, val):
    index = -1
    for i in xrange(0,len(array)):
        if val <= array[i]:
            index = i
            break
    if index == -1:
        array = np.append(array, val)
    else:
        array = np.insert(array, index, val)
    return array

def remove_particle_plot(file_name):
    ext = file_name.split(".txt")
    fos = open(ext[0]+'_removed.txt', 'r')
    lines = fos.readlines()
    eject = np.zeros(0)
    collide = np.zeros(0)
    time = np.zeros(0)
    Ntot=0
    Ncoll=0
    Nej=0
    for line in lines:
        split = line.split(',')
        val = float(split[1])
        if split[0] == "Collision":
            collide = sort_index(collide, val)
            Ncoll += 1
        elif split[0] == "Ejection":
            eject = sort_index(eject, val)
            Nej += 1
        time = sort_index(time, val)
        Ntot += 1
    axes[3].plot(time, cdf(time)*Ntot, color='red', label='Avg. tot. removed particles')
    axes[3].plot(collide, cdf(collide)*Ncoll, color='green', linestyle='--',label='Avg. collided particles')
    axes[3].plot(eject, cdf(eject)*Nej, color='blue',linestyle='-.',label='Avg. ejected particles')
    axes[3].set_ylabel('removed particles')
    axes[3].set_xlabel('time')
    axes[3].legend(loc='upper left',prop={'size':10})
    axes[0].plot(collide, cdf(collide)*Ncoll*1e-12, color='brown', label='theory energy')

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

diagnostics = 1
plot_dr_or_com = 1  #0 = dr, 1 = com

fos = open(file_name, 'r')
data = np.loadtxt(fos, delimiter=',')

signed_energy = 0
reverse_xaxis = 0
removed_particles = 1
mom_or_com = 0  #plot momentum or com

ms=3
if diagnostics == 1:
    time,dE,Lang,Llin,com = data[:,0],data[:,1],data[:,9],data[:,10],data[:,11]
    fig, axes = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(10,10))
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
#momentum
    if mom_or_com == 1:
        if signed_energy == 1:
            a=[Lang>0]
            b=[Lang<0]
            c=[Llin<0]
            d=[Llin>0]
            axes[2].plot(time[a],abs(Lang[a]), 'o', ms=ms, markeredgecolor='none',label='Angular Momentum pos')
            axes[2].plot(time[b],abs(Lang[b]), 'o', ms=ms, markeredgecolor='none',label='Angular Momentum neg')
            axes[2].plot(time[c],abs(Llin[c]), 'o', ms=ms, markeredgecolor='none',label='Linear Momentum neg')
            axes[2].plot(time[d],abs(Llin[d]), 'o', ms=ms, markeredgecolor='none',label='Linear Momentum pos')
        else:
            axes[2].plot(time,abs(Lang), 'o', ms=ms, markeredgecolor='none',label='Angular Momentum')
            axes[2].plot(time,abs(Llin), 'o', ms=ms, markeredgecolor='none',label='Linear Momentum')
        axes[2].set_ylabel('momentum')
        axes[2].legend(loc='upper left',prop={'size':10})
        axes[2].set_yscale('log')
    else:
        axes[2].plot(time,com, 'o', ms=ms, markeredgecolor='none',label='')
        axes[2].plot(time,1e-11*time**1.5, color='black',label='t^1.5')
        axes[2].plot(time,1e-11*time**2, color='red',label='t^2')
        axes[2].set_xlabel('simulation time (yrs)')
        axes[2].set_ylabel('dcom')
        axes[2].set_yscale('log')
        axes[2].legend(loc='lower left',prop={'size':10})
    if removed_particles == 1:
        remove_particle_plot(file_name)

#HSR
    #axes[2].plot(time,data[:,10], 'o', ms=ms, markeredgecolor='none')
    #axes[2].set_ylabel('HSR')
#com

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
