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

def theory(system_params):
    a,Ms,Me,mp,dt,HSR = system_params   #system parameters
    rhHSR = HSR*a*(Me/(3*Ms))**(1./3.)
    a2 = a*a
    M3 = Me*mp*Ms
    tau2 = dt*dt / 12.
    #theory
    term1 = M3/(a*rhHSR**3)
    term2 = Me*Me*mp/(rhHSR**4)
    term3 = -M3/(2*rhHSR*((a2 - rhHSR**2)**1.5))    #minor term, neg^x returns invalue value
    termp = 3*M3/((a*rhHSR**3))
    theory = (term1 + term2 + term3)*tau2
    return theory

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

runs = ['HYBARIDcoll_Np500','HYBARIDcoll2xR_Np500','HYBARIDcoll4xR_Np500']
names = ['R','2R','4R']
colordark = ['darkgreen','darkblue','darkred']
colorlight = ['lightgreen','dodgerblue','salmon']
theoryparams = [0.5,1,1e-5,1e-8,0.001/(2*np.pi),6] #a,M*,ME,mp,dt,HSR

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
    axes1 = plt.subplot(gs[1],sharex=axes0)
    axes2 = plt.subplot(gs[2],sharex=axes0)
    axes = [axes0,axes1,axes2]
    plt.subplots_adjust(hspace = 0.3)

##############################################
for k in xrange(0,len(runs)):
    N_files = 0
    dirP = str(sys.argv[1])
    files = glob.glob(dirP+runs[k]+'_sd*.txt')
    N = len(files)
    i=0
    while i < N:
        f = files[i]
        string = f.split("_")
        if string[-1]=="removed.txt":
            files.remove(files[i])
            N -= 1
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
    Eavg = np.zeros(n_it)
    time = np.zeros(n_it)
    dcom = np.zeros(n_it)
    for i in xrange(0,n_it):
        vals_for_med = np.zeros(N_files)
        valsdcom = np.zeros(N_files)
        for j in range(0,N_files):
            split = data[j][i].split(",")
            vals_for_med[j] = float(split[1])      #energy
            #vals_for_med[j] = abs(float(split[9]))       #ang. mom.
            E[j][i] = vals_for_med[j]
            valsdcom[j] = float(split[11])
        Eavg[i] = np.median(vals_for_med)
        time[i] = float(split[0])
        dcom[i] = np.median(valsdcom)

    #plotting
    j=0
    if plot_only_avg == 0:
        for i in xrange(0,N_files):
            axes[0].plot(time,E[i], '.', color=colorlight[k], alpha=alpha)
    axes[0].plot(time, Eavg, '.', markeredgecolor='none', color=colordark[k], label=names[k])
    if POE == 0:
        axes[1].plot(time, dcom, '.', markeredgecolor='none', color=colordark[k])

    if plot_removed_particles == 1 and POE == 0: #removed particles
        dirP = str(sys.argv[1])
        files = glob.glob(dirP+runs[k]+'_sd*_removed.txt')
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
                val = float(split[1])
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

        #error growth theory
        scale_factor = 50000 #dE /(dE/E)
        HSR_err = theory(theoryparams)
        #HSR_err *= scale_factor

        #plotting
        if k==0:
            axes[2].plot(time, cdf(time)*Ntot, color=colordark[k], label='Avg. tot. removed particles')
            axes[2].plot(collide, cdf(collide)*Ncoll, colordark[k], linestyle='--',label='Avg. collided particles')
            axes[2].plot(eject, cdf(eject)*Nej, colordark[k],linestyle='-.',label='Avg. ejected particles')
            #axes[0].plot(collide, cdf(collide)*HSR_err, colordark[k],linestyle='-.',label='Error growth Theory')
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
#axes[0].set_ylabel('dL/L(0)', fontsize=fontsize)
axes[0].set_yscale('log')
axes[0].set_xscale('log')
axes[0].set_xlim([0.5,time[-1]])
if POE == 1:
    axes[0].set_xlabel('time (years)', fontsize=fontsize)
else:
    axes[1].plot(time,1e-13*time**1.5, color='black',label='t^1.5')
    axes[1].plot(time,1e-13*time**2, color='red',label='t^2')
    axes[1].set_ylabel('Averaged COM$_i$ - COM$_0$', fontsize=fontsize)
    axes[1].set_yscale('log')
    axes[1].legend(loc='upper left',prop={'size':10})
    axes[2].set_ylabel('removed particles', fontsize=fontsize)
    #axes[2].set_yscale('log')
    axes[2].set_xscale('log')
    axes[2].set_xlabel('time (years)', fontsize=fontsize)
    if plot_removed_particles == 1:
        axes[2].legend(loc='lower left',prop={'size':10})
print 'Preparing PDF'
plt.savefig(dirP+'energy_avgruns_'+outputname+'.png')
plt.show()