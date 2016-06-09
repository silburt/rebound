#This macro calculates the average energy for a set of rebound runs.
#dr = position difference between mini and global (for the massive bodies) just before the mini updates the global positions.
#removed = CDFs of removed particles (be it due to ejection or collision)

import glob
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pylab as pl
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import sys
import os

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

def cdf(array):
    return np.arange(1,len(array)+1)/float(len(array))

runs = ['Vanilla_powerlaw']
names = ['Kirsh']
colordark = ['darkgreen','darkblue','darkred']
colorlight = ['lightgreen','dodgerblue','salmon']

alpha_disk = [3,1,1]
a_ini = 25

outputname = runs[0]+'avg'

plot_removed_particles = 1
plot_only_avg = 0
POE = 0         #plot only energy

ms = 0.25
alpha = 0.4
fontsize = 14

if POE == 1:
    fig = plt.figure(figsize=(8, 6))
    gs = gridspec.GridSpec(1, 1)
    axes0 = plt.subplot(gs[0])
    axes = [axes0]
else:
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 1])
    axes0 = plt.subplot(gs[0])
    #axes1 = plt.subplot(gs[1],sharex=axes0)
    axes1 = plt.subplot(gs[1])
    axes2 = plt.subplot(gs[2],sharex=axes0)
    axes = [axes0,axes1,axes2]
    plt.subplots_adjust(hspace = 0.3)

##############################################
for k in xrange(0,len(runs)):
    N_files = 0
    dirP = str(sys.argv[1])
    files = glob.glob(dirP+runs[k]+'*.txt')
    N = len(files)
    i=0
    while i < N:
        f = files[i]
        string = f.split("_")
        if string[-1]=="removed.txt" or string[-1]=="elapsedtime.txt":
            files.remove(files[i])
            N -= 1
        else:
            i += 1
    data = []
    n_it = 10e10
    for f in files:
        try:
            ff = open(f, 'r')
            lines = ff.readlines()
            length = len(lines)
            if length < n_it:   #need to find array with shortest length
                n_it = length
            data.append(lines)
            N_files += 1
        except:
            print 'couldnt read in data file '+f

    E = np.zeros(shape=(N_files,n_it))
    a = np.zeros(shape=(N_files,n_it))
    aavg = np.zeros(n_it)
    Eavg = np.zeros(n_it)
    time = np.zeros(n_it)
    for i in xrange(0,n_it):
        vals_for_medE = np.zeros(N_files)
        vals_for_meda = np.zeros(N_files)
        for j in range(0,N_files):
            split = data[j][i].split(",")
            vals_for_medE[j] = float(split[1])
            vals_for_meda[j] = float(split[2])
            E[j][i] = vals_for_medE[j]
            a[j][i] = vals_for_meda[j]
        Eavg[i] = np.median(vals_for_medE)
        aavg[i] = np.median(vals_for_meda)
        time[i] = float(split[0])/6.283185

    #plotting
    j=0
    if plot_only_avg == 0:
        for i in xrange(0,N_files):
            axes[0].plot(time,abs(E[i]), '.', color=colorlight[k], alpha=alpha)
            axes[1].plot(time,a[i], '.', color=colorlight[k], alpha=alpha)
    axes[0].plot(time, abs(Eavg), '.', markeredgecolor='none', color=colordark[k], label=names[k])
    axes[1].plot(time, aavg, '.', markeredgecolor='none', color=colordark[k], label=names[k])

    #semi-major axis migration theory - Abstract of Kirsh
    P = np.sqrt(4* np.pi**2 * a_ini**3)
    dt = time[1:] - time[0:-1]
    #a = aavg[1:] - dt*factor*aavg[1:]**(3 - alpha_disk[k])
    #axes[1].plot(time[0:-1], a, 'k.',label='theory (alpha='+str(alpha_disk[k])+')')
    dadta = 4*np.pi*(1.913e-4)*a_ini**(3 - alpha_disk[k])/(P*a_ini)
    #a = a_ini*(1 - (alpha_disk[k] + 0.5)*dadta*time)**(-1./(alpha_disk[k] + 0.5))
    factor = (1.913e-4)*4*(np.pi/P)
    a = a_ini - factor*time*a_ini**(3-alpha_disk[k])/25
    #a = (factor*time*(2 - alpha_disk[k]) + a_ini**(alpha_disk[k] - 2))**(1./(alpha_disk[k]-2))
    print a
    axes[1].plot(time, a, 'b.',label='theory (alpha='+str(alpha_disk[k])+')')

    if plot_removed_particles == 1 and POE == 0: #removed particles
        dirP = str(sys.argv[1])
        files = glob.glob(dirP+runs[k]+'*_removed.txt')
        eject = np.zeros(0)
        collide = np.zeros(0)
        time = np.zeros(0)
        Ntot=0
        Ncoll=0
        Nej=0
        for f in files:
            ff = open(f, 'r')
            lines = ff.readlines()
            for line in lines:
                split = line.split(',')
                val = float(split[1])/6.283185
                if split[0] == "Collision":
                    collide = sort_index(collide, val)
                    Ncoll += 1
                elif split[0] == "Ejection":
                    eject = sort_index(eject, val)
                    Nej += 1
                time = sort_index(time, val)
                Ntot += 1

        Ntot /= len(files)
        Ncoll /= len(files)
        Nej /= len(files)

        #plotting
        if k==0:
            axes[2].plot(time, cdf(time)*Ntot, color=colordark[k], label='Avg. tot. removed particles')
            axes[2].plot(collide, cdf(collide)*Ncoll, colordark[k], linestyle='--',label='Avg. collided particles')
            axes[2].plot(eject, cdf(eject)*Nej, colordark[k],linestyle='-.',label='Avg. ejected particles')
        else:
            axes[2].plot(time, cdf(time)*Ntot, color=colordark[k])
            axes[2].plot(collide, cdf(collide)*Ncoll, colordark[k], linestyle='--')
            axes[2].plot(eject, cdf(eject)*Nej, colordark[k],linestyle='-.')

##############################################
#Final plotting stuff
axes[0].plot(time,0.8e-11*time**(0.5),color='black', label='t$^{ 0.5}$')
#axes[0].plot(time,2e-14*time,color='red', label='t')
axes[0].legend(loc='upper left',prop={'size':10}, numpoints=1, markerscale=3)
axes[0].set_ylabel('dE/E(0)', fontsize=fontsize)
axes[0].set_yscale('log')
axes[0].set_xscale('log')
axes[0].set_xlim([20,time[-1]])
if POE == 1:
    axes[0].set_xlabel('time (years)', fontsize=fontsize)
else:
    axes[1].set_ylabel('a(t)', fontsize=fontsize)
    axes[1].legend(loc='upper right',prop={'size':10})
    axes[1].set_ylim([3,7])
    axes[2].set_ylabel('removed particles', fontsize=fontsize)
    axes[2].set_yscale('log')
    axes[2].set_xscale('log')
    axes[2].set_xlabel('time (years)', fontsize=fontsize)
    if plot_removed_particles == 1:
        axes[2].legend(loc='lower left',prop={'size':10})
print 'Preparing PDF'
plt.savefig(dirP+'energy_avgruns_'+outputname+'.png')
plt.show()