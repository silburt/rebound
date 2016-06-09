#This macro shows the evolution of semimajor axis, eccentricity and inclination for different snapshots in time. Done in a serial fashion.

import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import re

ext = '_ei0.txt'

dir = sys.argv[1]
runs = glob.glob(dir+'*'+ext)
Nruns = len(runs)

for i in xrange(0,Nruns):
    snapshot_name = re.sub('0\.txt$', '', runs[i])
    snapshots = glob.glob(snapshot_name+'*.txt')
    Nsnapshots = len(snapshots)
    if Nsnapshots > 1:
        fig, axes = plt.subplots(nrows=Nsnapshots, ncols=1, figsize=(12,12), sharex=True)
        for j in xrange(0,Nsnapshots):
            fos = open(snapshots[j], 'r')
            time, id, a, e, inc, rdist, m = np.loadtxt(fos, delimiter=',',unpack=True)
            #ei = (e*e + 2*inc*inc)**(0.5)
            axes[j].plot(a,e,'.')
            axes[j].plot(a[0],e[0],'o',color='red')
            axes[j].plot(a[1],e[1],'o',color='orange')
            axes[j].set_title('t = '+str(time[0]/6.283)+' years, N='+str(len(a)))
            axes[j].yaxis.set_ticks(np.arange(0,max(e),0.1))
        #axes[Nsnapshots-1].set_ylabel('(e$^2$ + 2I$^2$)$^{0.5}$')
        axes[Nsnapshots-1].set_ylabel('$e$')
        axes[Nsnapshots-1].set_xlabel('(a (AU)')
        plt.savefig(snapshot_name+'.png')
        print 'done file',i
