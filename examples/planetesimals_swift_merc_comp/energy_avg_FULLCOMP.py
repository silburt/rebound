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
from subprocess import call

#assuming M=G=1!!
def calc_a(x,y,z,vx,vy,vz):
    v2 = (vx**2 + vy**2 + vz**2)
    d = (x**2 + y**2 + z**2)**0.5
    return -1./(v2 - 2./d)

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
            lines = lines[9:-1]     #for mercury
            length = len(lines)
            if length < n_it:   #need to find array with shortest length
                n_it = length
            data.append(lines)
            N_files += 1
        except:
            print 'couldnt read in data file '+f
    return data, n_it, N_files

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
            elapsed.append(float(lines[0].split()[-1])/3600.)
    return elapsed

def cdf(array):
    return np.arange(1,len(array)+1)/float(len(array))

swifter = 0
mercury = 0

fig = plt.figure(figsize=(10, 10))
gs = gridspec.GridSpec(4, 1, height_ratios=[2, 1, 1, 1])
axes = [plt.subplot(gs[0]),plt.subplot(gs[1]),plt.subplot(gs[2]),plt.subplot(gs[3])]
plt.subplots_adjust(hspace = 0.3)
ms = 0.25
alpha = 0.4
counter = 1

##############################################
#HERMES
dirP = str(sys.argv[1])
#name = ['t5e7','Earth','ias15_', 'whfast']
#name = ['t1e+06']
name = ['t5e+07']
outname = ['HERMES:Neptune','HERMES:Earth-sized','IAS15:Neptune', 'WHFAST:Neptune']
color_back = ['lightgreen','violet','yellow','navajowhite']
color_main = ['darkgreen','darkviolet','olive','darkorange']
ETarr = []
for ii in range(len(name)):
    files = glob.glob(dirP+name[ii]+'*.txt')
    N = len(files)
    i=0
    while i < N:    #just want the main .txt files
        f = files[i]
        string = f.split("_")
        if ("removed" in string[-1]) or ("elapsed" in string[-1]) or ("ET" in string[-1]):
            files.remove(files[i])
            N -= 1
        else:
            i += 1

    Energy = []
    Npar = []
    aa = []
    for f in files:
        #time, E, Np, miniN, mini_active, elapsedtime, a, HSF = np.loadtxt(f,delimiter=',',unpack=True)
        #time, E, Np, miniN, mini_active, a, HSF = np.loadtxt(f,delimiter=',',unpack=True)
        time, E, Np, miniN, mini_active, a, HSF, stepsminiactive = np.loadtxt(f,delimiter=',',unpack=True)
        axes[0].plot(time/(2*np.pi),E, '.', color=color_back[ii], alpha=alpha)
        aa.append(a)
        #axes[2].plot(time/(2*np.pi),a,'.',color=color_back[ii],alpha=alpha)
        Energy.append(E)
        Npar.append(Np)

    try:
        Energy = np.asarray(zip(*Energy))
        Npar = np.asarray(zip(*Npar))
        aa = np.asarray(zip(*aa))
        axes[0].plot(time/(2*np.pi),np.mean(Energy,axis=1), '.', markeredgecolor='none', color=color_main[ii], label=outname[ii])
        axes[0].plot(time/(2*np.pi), 1e-10*time**0.5, color='black')
        axes[0].plot(time/(2*np.pi), 0.4e-16*time**0.5, color='black')
        axes[2].plot(time/(2*np.pi),np.mean(aa,axis=1), '.', markeredgecolor='none', color=color_main[ii])
        axes[3].plot(time/(2*np.pi),np.mean(Npar,axis=1), '.', markeredgecolor='none', color=color_main[ii])
    except:
        axes[0].plot(time/(2*np.pi),E, '.', color=color_back[ii], alpha=alpha,label=outname[ii])
        axes[2].plot(time/(2*np.pi),a,'.',color=color_back[ii],alpha=alpha,label=outname[ii])
        print "each file has different dimensions, probably IAS15"
    
    #Elapsed time
    files = glob.glob(dirP+name[ii]+'*_elapsedtime.txt')
    times_H = get_times(files,'','H')
    axes[1].plot(times_H, np.ones(len(times_H))*counter, 'o', markeredgecolor='none', ms=10, color=color_main[ii])
    counter += 1
    ETarr.append(outname[ii])

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
    data, n_it, N_files = get_data(files,'/eo.txt')
    E = np.zeros(shape=(N_files,n_it))
    Eavg = np.zeros(n_it)
    time = np.zeros(n_it)
    vals_for_med = np.zeros(N_files)
    nump = np.zeros(n_it)
    vals_for_nump = np.zeros(N_files)

    print 'calculating avg. energy, new v.'
    for i in xrange(0,n_it):
        split = data[0][i].split()
        vals_for_med[0] = float(split[4])
        E[0][i] = vals_for_med[0]
        t=float(split[1])
        for j in range(1,N_files):
            split = data[j][i].split()
            if float(split[1]) == t:
                vals_for_med[j] = float(split[4])
                #vals_for_nump[j] = int(split[7])
                E[j][i] = vals_for_med[j]
            else:
                print 'problem, times dont match'
        Eavg[i] = np.median(vals_for_med)
        time[i] = float(split[1])
        #nump[i] = np.median(vals_for_nump)

    for i in xrange(0,N_files):
        axes[0].plot(time,E[i], '.', color='salmon', alpha=alpha)
    axes[0].plot(time, Eavg, '.', markeredgecolor='none', color='darkred', label='MERCURY Avg.')
    #axes[3].plot(time, nump, '.', markeredgecolor='none', color='darkred')
'''
    #get median semi-major axis
    r = []
    day2yr2pi = 365./(2.*np.pi)
    for i in range(len(files)):
        f = files[i]
        try:
            time2, mass, x, y, z, vx, vy, vz = np.loadtxt(f+'/BODY1.aei',skiprows=4,unpack=True)
        except:
            call('(cd %s && gfortran -o element6 element6.for && ./element6)'%f,shell=True)
            time2, mass, x, y, z, vx, vy, vz = np.loadtxt(f+'/BODY1.aei',skiprows=4,unpack=True)
        r.append(calc_a(x,y,z,vx*day2yr2pi,vy*day2yr2pi,vz*day2yr2pi))
    r = np.asarray(zip(*r))
    axes[2].plot(time2,np.mean(r,axis=1),'.', color='salmon',alpha=alpha)

    #Elapsed time
    times_M = get_times(files,'/ET.txt','M')
    axes[1].plot(times_M, np.ones(len(times_M))*counter, 'o',markeredgecolor='none', ms=10,color='salmon')
    counter += 1
    ETarr.append('MERCURY')
'''
##############################################
#Final plotting stuff
fontsize=10
xmin = 1
axes[0].legend(loc='upper left',prop={'size':7}, numpoints=1, markerscale=3)
axes[0].set_ylabel('Fractional Energy Error', fontsize=fontsize)
axes[0].set_xlabel('Time (Years)', fontsize=fontsize)
axes[0].set_yscale('log')
axes[0].set_ylim([1e-16, 1])
axes[0].set_xscale('log')
axes[0].set_xlim([xmin,time[-1]])
axes[1].set_ylabel('Elapsed time (hours)')
axes[1].set_xlabel('Elapsed Simulation Time (Hours)',fontsize=fontsize)
axes[1].set_yticks(range(1,7))
axes[1].set_yticklabels(ETarr,fontsize=fontsize)
axes[1].set_ylim([0.5,6.5])
axes[2].set_ylabel('Semi-major axis')
axes[2].set_xlabel('Time (Years)', fontsize=fontsize)
axes[3].set_ylabel('Number of particles')
axes[3].set_xscale('log')
axes[3].set_xlabel('Time (Years)', fontsize=fontsize)
axes[3].set_xlim([xmin,time[-1]])

print 'Preparing PDF'
#plt.savefig(pp, format='pdf')
plt.savefig(dirP+'energy_avg_FULL.png')
#plt.show()

#hermes old
'''
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
    axes[0].plot(time,E[i], '.', color=color_back[ii], alpha=alpha)
    axes[0].plot(time, Eavg, '.', markeredgecolor='none', color=color_main[ii], label=outname[ii])
    axes[0].plot(time, 1e-10*time**0.5, color='black')
    axes[0].plot(time, 0.4e-16*time**0.5, color='black')
    '''

#mercury old
'''
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
    '''
