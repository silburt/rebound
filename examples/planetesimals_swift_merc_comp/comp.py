#This is the revamped plot going in the HERMES paper.

import glob
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pylab as pl
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import sys
import os
from subprocess import call

def get_times(files, ext, integ):
    elapsed = []
    for f in files:
        ff = open(f+ext, 'r')
        lines = ff.readlines()
        if integ == 'H':
            elapsed.append(float(lines[1].split()[-2])/3600.)
        elif integ == 'S':
            elapsed.append((float(lines[1].split()[-1]) - float(lines[0].split()[-1]))/3600.)
        elif integ == 'M':
            elapsed.append(float(lines[0].split()[-2])/3600.)
    return elapsed

swifter = 0
mercury = 1

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
axes = [plt.subplot(gs[0]),plt.subplot(gs[1])]
plt.subplots_adjust(hspace = 0.3)
ms = 0.25
alpha = 0.4
counter = 1
ETarr = []

##############################################
#SWIFTER
if swifter == 1:
    dir = '../../../swifter/example/input_files/'
    files = [x[0] for x in os.walk(dir)][1:]
    dE_arrS, t_arrS = [], []
    for f in files:
        time, dE, N, offset = np.loadtxt(open(f+'/energyoutput.txt','r'), unpack=True)
        t_arrS.append(time), dE_arrS.append(dE)
        axes[0].plot(time,dE, '.', color='dodgerblue', alpha=alpha)
    dE_meanS = np.median(np.asarray(zip(*dE_arrS)),axis=1)
    t_meanS = np.median(np.asarray(zip(*t_arrS)),axis=1)
    axes[0].plot(t_meanS, dE_meanS, '.', markeredgecolor='none', color='darkblue', label='SyMBA Avg.')

    #Elapsed time
    times_S = get_times(files,'/elapsed_time.txt','S')
    axes[1].plot(times_S, np.ones(len(times_S))*counter, 'o', markeredgecolor='none', ms=10, color='dodgerblue')
    counter += 1
    ETarr.append('SyMBA')

##############################################
#MERCURY
if mercury == 1:
    print '...Finished Swifter, doing Mercury now...'
    dir = '../../../mercury6/input_files/'
    files = [x[0] for x in os.walk(dir)][1:]
    dE_arrM, t_arrM = [], []
    for f in files:
        data = np.genfromtxt(f+'/eo.txt', delimiter=None, dtype=None, skip_header=2, skip_footer=2)
        junk,time,junk,junk,dE,junk,dL,N = zip(*data.T)
        #time *= 0.0172142
        t_arrM.append(time), dE_arrM.append(dE)
        axes[0].plot(time,dE, '.', color='salmon', alpha=alpha)

    dE_meanM = np.median(np.asarray(zip(*dE_arrM)),axis=1)
    t_meanM = np.median(np.asarray(zip(*t_arrM)),axis=1)
    axes[0].plot(t_meanM, dE_meanM, '.', markeredgecolor='none', color='darkred', label='MERCURY Avg.')

    #Elapsed Time
    times_M = get_times(files,'/ET.txt','M')
    axes[1].plot(times_M, np.ones(len(times_M))*counter, 'o',markeredgecolor='none', ms=10,color='salmon')
    counter += 1
    ETarr.append('MERCURY')

##############################################
#HERMES
print '...Doing Hermes now...'
dirH = str(sys.argv[1])
files = glob.glob('%s*_elapsedtime.txt'%dirH)
dE_arr, t_arr = [], []
for f in files:
    f2 = f.split('_elapsedtime.txt')[0]+'.txt'
    time, dE, N, miniN, miniactive, a, HSF, MAS = np.loadtxt(open(f2,'r'), delimiter=',', unpack=True)
    time /= 2*np.pi
    t_arr.append(time), dE_arr.append(dE)
    axes[0].plot(time,dE, '.', color='lightgreen', alpha=alpha)

dE_mean = np.median(np.asarray(zip(*dE_arr)),axis=1)
t_mean = np.median(np.asarray(zip(*t_arr)),axis=1)
times_H = get_times(files,'','H')
axes[0].plot(t_mean, dE_mean, '.', markeredgecolor='none', color='darkgreen', label='HERMES Avg.')
axes[0].plot(t_mean, 5e-10*t_mean**0.5, '.', markeredgecolor='none', color='black', label='sqrt(t)')
axes[1].plot(times_H, np.ones(len(times_H))*counter, 'o', markeredgecolor='none', ms=10, color='lightgreen')
counter += 1
ETarr.append('HERMES')

##############################################
#Final plotting stuff
fontsize=10
xmin = 1
axes[0].legend(loc='upper left',prop={'size':7}, numpoints=1, markerscale=3)
axes[0].set_ylabel('Fractional Energy Error', fontsize=fontsize)
axes[0].set_xlabel('Time (Years)', fontsize=fontsize)
axes[0].set_yscale('log')
axes[0].set_ylim([1e-12, 1e-3])
axes[0].set_xscale('log')
axes[0].set_xlim([xmin,time[-1]])
axes[1].set_ylabel('Elapsed time (hours)')
axes[1].set_xlabel('Elapsed Simulation Time (Hours)',fontsize=fontsize)
axes[1].set_yticks(range(1,7))
axes[1].set_yticklabels(ETarr,fontsize=fontsize)
axes[1].set_ylim([0.5,6.5])

print 'Preparing PDF'
#plt.savefig(pp, format='pdf')
plt.savefig(dirH+'comp.png')
#plt.show()
