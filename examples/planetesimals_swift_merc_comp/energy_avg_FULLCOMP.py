#This macro calculates the average energy for a set of PLANETESIMAL, SWIFTER and MERCURY runs and plots them all on one sexy figure.
#This plot is going in the paper

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

def get_data(files, ext):
    data = []
    n_it = 10e10
    N_files = 0
    for f in files:
        try:
            ff = open(f+ext, 'r')
            lines = ff.readlines()
            length = len(lines)
            if length < n_it:   #need to find array with shortest length
                n_it = length
            data.append(lines)
            N_files += 1
        except:
            print 'couldnt read in data file '+f
    return data, n_it, N_files

def get_times(files, ext):
    elapsed = []
    for f in files:
        ff = open(f+ext, 'r')
        lines = ff.readlines()
        if ext == '':
            elapsed.append(float(lines[1].split()[-2])/3600.)
        else:
            elapsed.append((float(lines[1].split()[-1]) - float(lines[0].split()[-1]))/3600.)
    return elapsed

def cdf(array):
    return np.arange(1,len(array)+1)/float(len(array))

swifter = 1
mercury = 1

fig = plt.figure(figsize=(8, 8))
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
axes0 = plt.subplot(gs[0])
axes1 = plt.subplot(gs[1])
axes = [axes0,axes1]
plt.subplots_adjust(hspace = 0.1)
ms = 0.25
alpha = 0.4

##############################################
#PLANETESIMAL
dirP = str(sys.argv[1])
files = glob.glob(dirP+'*.txt')
N = len(files)
i=0
while i < N:    #just want the main .txt files
    f = files[i]
    string = f.split("_")
    if ("removed" in string[-1]) or ("elapsed" in string[-1]):
        files.remove(files[i])
        N -= 1
    else:
        i += 1

data, n_it, N_files = get_data(files,'')
E = np.zeros(shape=(N_files,n_it))
Eavg = np.zeros(n_it)
time = np.zeros(n_it)
vals_for_med = np.zeros(N_files)
for i in xrange(0,n_it):
    for j in range(0,N_files):
        split = data[j][i].split(",")
        vals_for_med[j] = float(split[1])
        E[j][i] = vals_for_med[j]
    Eavg[i] = np.median(vals_for_med)
    time[i] = float(split[0])

for i in xrange(0,N_files):
    axes[0].plot(time,E[i], '.', color='lightgreen', alpha=alpha)
axes[0].plot(time, Eavg, '.', markeredgecolor='none', color='darkgreen', label='HERMES avg.')

#Elapsed time
files = glob.glob(dirP+'*_elapsedtime.txt')
times_H = get_times(files,'')

##############################################
#SWIFTER
if swifter == 1:
    print '...Finished Planetesimal, doing Swifter now...'
    dir = '../../../swifter/example/input_files/'
    files = [x[0] for x in os.walk(dir)][1:]
    data, n_it, N_files = get_data(files,'/energyoutput.txt')
    E = np.zeros(shape=(N_files,n_it))
    Eavg = np.zeros(n_it)
    time = np.zeros(n_it)
    vals_for_med = np.zeros(N_files)
    print 'calculating avg energy'
    for i in xrange(1,n_it):
        for j in range(0,N_files):
            split = data[j][i].split()
            vals_for_med[j] = float(split[1])
            E[j][i] = vals_for_med[j]
        Eavg[i] = np.median(vals_for_med)
        time[i] = float(split[0])

    j=0
    for i in xrange(0,N_files):
        axes[0].plot(time,E[i], '.', color='dodgerblue', alpha=alpha)
    axes[0].plot(time, Eavg, '.', markeredgecolor='none', color='darkblue', label='SyMBA Avg.')

#Elapsed time
times_S = get_times(files,'/elapsed_time.txt')

##############################################
#MERCURY
if mercury == 1:
    print '...Finished Swifter, doing Mercury now...'
    dir = '../../../mercury6/input_files/'
    files = [x[0] for x in os.walk(dir)][1:]
    data, n_it, N_files = get_data(files,'/eo.txt')
    E = np.zeros(shape=(N_files,n_it))
    Eavg = np.zeros(n_it)
    time = np.zeros(n_it)
    vals_for_med = np.zeros(N_files)

    print 'calculating avg energy'
    for i in xrange(0,n_it):
        split = data[0][i].split()
        vals_for_med[0] = float(split[1])
        E[0][i] = vals_for_med[0]
        t=float(split[0])
        for j in range(1,N_files):
            split = data[j][i].split()
            if float(split[0]) == t:
                vals_for_med[j] = float(split[1])
                E[j][i] = vals_for_med[j]
            else:
                print 'problem, times dont match'
        Eavg[i] = np.median(vals_for_med)
        time[i] = float(split[0])*0.01721420632

    for i in xrange(0,N_files):
        axes[0].plot(time,E[i], '.', color='salmon', alpha=alpha)
    axes[0].plot(time, Eavg, '.', markeredgecolor='none', color='darkred', label='MERCURY Avg.')

#Elapsed time
times_M = get_times(files,'/elapsed_time.txt')

##############################################
#Final plotting stuff
axes[0].legend(loc='upper left',prop={'size':10}, numpoints=1, markerscale=3)
axes[0].set_ylabel('dE/E(0)', fontsize=16)
axes[0].set_yscale('log')
axes[0].set_ylim([1e-12, 1e-3])
axes[0].set_xscale('log')
axes[0].set_xlim([0.5,time[-1]])

axes[1].plot(times_H, np.ones(len(times_H)), 'o', markeredgecolor='none', color='lightgreen')
axes[1].plot(times_S, np.ones(len(times_S))*2, 'o', markeredgecolor='none', color='dodgerblue')
axes[1].plot(times_M, np.ones(len(times_M))*3, 'o',markeredgecolor='none', color='salmon')
axes[1].set_xlabel("Elapsed Simulation Time (hours)")
axes[1].set_yticks(range(1,4))
axes[1].set_yticklabels(['HERMES','SyMBA','MERCURY'])
axes[1].set_ylim([0.5,3.5])
print 'Preparing PDF'

#from matplotlib.backends.backend_pdf import PdfPages
#pp = PdfPages(dirP+'energy_avg_FULL.pdf')
#plt.savefig(pp, format='pdf')
plt.savefig(dirP+'energy_avg_FULL.pdf', format='pdf')
plt.show()




#used to be for HERMES
'''
    #Removed Particles
    dirP = str(sys.argv[1])
    files = glob.glob(dirP+'*_removed.txt')
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
    axes[1].plot(time, cdf(time)*Ntot, 'k', label='Avg. tot. removed particles (N='+str(Ntot)+')')
    axes[1].plot(collide, cdf(collide)*Ncoll, 'k--',label='Avg. collided particles (N='+str(Ncoll)+')')
    axes[1].plot(eject, cdf(eject)*Nej, 'k-.',label='Avg. ejected particles (N='+str(Nej)+')')
    axes[1].legend(loc='upper left',prop={'size':10})
    axes[1].set_xlabel('time (years)', fontsize=16)
    axes[1].set_ylabel('CDF of removed particles', fontsize=16)
    '''
