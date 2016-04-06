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

outputname='CEplot_Np500_nocomnocollforceinterp'
runs = ['Np500_nocollnocomforceinterp']
#runs = ['Np500_forceinterpnocoll']
runmasses = [1e-8,1e-8]
names = ['Np500, mp=1e-8','Np5000, mp=1e-9']
colordark = ['darkgreen','darkblue','darkred']
colorlight = ['lightgreen','dodgerblue','salmon']

plot_CE = 0
plot_removed_particles = 0
plot_only_avg = 0

ms = 0.25
alpha = 0.4
fontsize = 14

fig = plt.figure(figsize=(10, 10))
if plot_CE == 1:
    gs = gridspec.GridSpec(4, 1, height_ratios=[1.5, 1, 1, 1])
    axes0 = plt.subplot(gs[0])
    axes1 = plt.subplot(gs[1],sharex=axes0)
    axes2 = plt.subplot(gs[2],sharex=axes0)
    axes4 = plt.subplot(gs[3],sharex=axes0)
    axes = [axes0,axes1,axes2,axes4]
else:
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
    files = glob.glob(dirP+'t*_'+runs[k]+'_sd*.txt')
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
    com_gp = np.zeros(n_it) #com global partial (i.e. just massive bodies)
    com_gf = np.zeros(n_it) #com global full (i.e. all particles)
    com_mp = np.zeros(n_it)
    com_mf = np.zeros(n_it)
    if plot_CE == 1:
        CE = np.zeros(n_it)
        CEmass = np.zeros(n_it)
    for i in xrange(0,n_it):
        vals_for_med = np.zeros(N_files)
        valscom_gp = np.zeros(N_files)
        valscom_gf = np.zeros(N_files)
        valscom_mp = np.zeros(N_files)
        valscom_mf = np.zeros(N_files)
        for j in range(0,N_files):
            split = data[j][i].split(",")
            vals_for_med[j] = float(split[1])
            E[j][i] = vals_for_med[j]
            valscom_gp[j] = float(split[9])
            valscom_gf[j] = float(split[10])
            valscom_mp[j] = float(split[11])
            valscom_mf[j] = float(split[12])
        Eavg[i] = np.median(vals_for_med)
        time[i] = float(split[0])
        com_gp[i] = np.median(valscom_gp)
        com_gf[i] = np.median(valscom_gf)
        com_mp[i] = np.median(valscom_mp)
        com_mf[i] = np.median(valscom_mf)
        if plot_CE == 1:
            CE[i] = int(split[5])
            CEmass[i] = CE[i]*runmasses[k]

    #plotting
    j=0
    if plot_only_avg == 0:
        for i in xrange(0,N_files):
            axes[0].plot(time,E[i], '.', color=colorlight[k], alpha=alpha)
    axes[0].plot(time, Eavg, '.', markeredgecolor='none', color=colordark[k], label=names[k])
    if k==0:
        com_gp = abs((com_gp - com_gp[0]))
        com_gf = abs((com_gf - com_gf[0]))
        com_mp = abs((com_mp - com_mp[0]))
        com_mf = abs((com_mf - com_mf[0]))
        axes[2].plot(time, com_gp, color=colordark[k], label=names[k]+' COM global massive bodies')
        axes[2].plot(time, com_gf, linestyle='--', color=colordark[k], label=names[k]+' COM global all bodies')
        axes[2].plot(time, com_mp, linestyle='-.', color=colordark[k], label=names[k]+' COM mini massive bodies')
        axes[2].plot(time, com_mf, linestyle=':', color=colordark[k], label=names[k]+' COM mini all bodies')
        axes[2].set_ylabel('COM$_i$ - COM$_0$', fontsize=fontsize)
        if plot_CE == 1:
            axes[3].plot(time, CE, color=colordark[k], label=names[k]+' CE')
            axes[3].plot(time, CEmass, color=colordark[k], linestyle='--', label=names[k]+' CE*m_planetesimal')
    else:
        axes[2].plot(time, com_gp, color=colordark[k])
        axes[2].plot(time, com_gf, linestyle='--', color=colordark[k])
        axes[2].plot(time, com_mp, linestyle='-.', color=colordark[k])
        axes[2].plot(time, com_mf, linestyle=':', color=colordark[k])
        if plot_CE == 1:
            axes[3].plot(time, CE, color=colordark[k])
            axes[3].plot(time, CEmass, color=colordark[k],linestyle='--')

    if plot_removed_particles == 1: #removed particles
        dirP = str(sys.argv[1])
        files = glob.glob(dirP+'t*_'+runs[k]+'_sd*_removed.txt')
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

        #plotting
        if k==0:
            axes[1].plot(time, cdf(time)*Ntot, color=colordark[k], label='Avg. tot. removed particles')
            axes[1].plot(collide, cdf(collide)*Ncoll, colordark[k], linestyle='--',label='Avg. collided particles')
            axes[1].plot(eject, cdf(eject)*Nej, colordark[k],linestyle='-.',label='Avg. ejected particles')
        else:
            axes[1].plot(time, cdf(time)*Ntot, color=colordark[k])
            axes[1].plot(collide, cdf(collide)*Ncoll, colordark[k], linestyle='--')
            axes[1].plot(eject, cdf(eject)*Nej, colordark[k],linestyle='-.')

##############################################
#Final plotting stuff
axes[0].plot(time,1e-11*time**(0.5),color='black', label='t$^{ 0.5}$')
axes[0].plot(time,1e-12*time,color='red', label='t')
axes[0].legend(loc='lower left',prop={'size':10}, numpoints=1, markerscale=3)
axes[1].legend(loc='upper left',prop={'size':10})
axes[1].set_ylabel('removed particles', fontsize=fontsize)
axes[0].set_ylabel('dE/E(0)', fontsize=fontsize)
axes[0].set_yscale('log')
axes[0].set_xscale('log')
axes[1].set_yscale('log')
axes[0].set_xlim([0.5,time[-1]])
axes[1].set_xscale('log')
if plot_removed_particles == 1:
    axes[2].legend(loc='lower right',prop={'size':10})
if plot_CE == 1:
    axes[3].set_yscale('log')
    axes[3].legend(loc='upper left',prop={'size':10})
    axes[3].set_ylabel('Close encounters', fontsize=fontsize)
    axes[3].set_xlabel('time (years)', fontsize=fontsize)
else:
    axes[2].set_xlabel('time (years)', fontsize=fontsize)
print 'Preparing PDF'
plt.savefig(dirP+'energy_avgruns_'+outputname+'.png')
plt.show()