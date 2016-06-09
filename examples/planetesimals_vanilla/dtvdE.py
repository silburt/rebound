#This macro makes a plot of dt vs. energy (and also dt vs. time)

import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import sys
import os
import re
twopi = 6.283185307

def sort(x, y):
    length = len(x)
    for i in xrange(0,length):
        for j in xrange(i,length):
            if x[j] < x[i]:
                tempx = x[j]
                tempy0 = y[j]
                x[j] = x[i]
                x[i] = tempx
                y[j] = y[i]
                y[i] = tempy0
    return x, y

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

#plot dt vs.
y_choice = int(sys.argv[2])      #0 = plot elapsed time, 1 = energy, 2=semi-major axis
x_choice = 2                    #0 = dt, 1 = HSR, 2 = powerlaw

dirP = sys.argv[1]
Navg = 20                       #number of points at end of .txt file to average energy over
xvals = []

#if x_choice == 0:   #dt
#    ext = 'dt'
#    xname = 'timestep (years)'
#    yvals = ['Np0']
#    path = '_th*_elapsedtime.txt'
#    labels = yvals
#    marker = '.'
#    title = 'dt: Integrating 2p, 2pl system for 1000 orbits'
#    HSR = 6
#if x_choice == 1:   #HSR
#    ext = 'HSR'
#    xname = 'HSR (hill radii)'
#    #yvals = ['Np500','Np5000']
#    yvals = ['Np0']
#    path = '_th*_elapsedtime.txt'
#    labels = yvals
#    marker = '.'
#    title = 'HSR: Integrating 2p, 2pl system for 1M orbits'
#    dt = np.array(1e-3)/twopi  #must be in yr/2pi format
if x_choice == 2:   #Np
    ext = ['powerlaw']
    xname = 'powerlaw of planetesimals'
    yvals = ['Np=5000, HSR=6, Mdisk = 10mp']
    path = '_Np5000_HSR6_Mpl10mp_elapsedtime.txt'
    labels = yvals
    marker = '.'
    title = ''
    HSR = 6
    dt = np.array(1e-3)/twopi #must be in yr/2pi format

a = np.array(1)             #AU
Ms = np.array(1)            #Solar mass
mp = np.array(1e-8)         #Solar mass
Me = np.array(5e-5)         #Solar mass
rh = a*(Me/(3*Ms))**(1./3.) #AU

leny = len(ext)
for i in xrange(0,leny):
    files = glob.glob(dirP+'Vanilla_'+ext[i]+'*'+path)
    files = sorted(files, key = natural_key)    #assuming files are sorted in correct order now
    lenf = len(files)
    ET = np.zeros(lenf)
    dE = np.zeros(lenf)
    a = np.zeros(lenf)
    xvals = np.zeros(lenf)
    for j in xrange(0,lenf):
        f = files[j]
        header = f.split("_")
        split1 = header[-5]
        split2 = split1.split(ext[i])
        xvals[j] = float(split2[1])
        ff = open(f, 'r')
        lines = ff.readlines()
        elapsed = lines[1]
        elapsed = elapsed.split()
        txtfile = f.split("_elapsedtime.txt")
        fff = open(txtfile[0]+".txt","r")
        lines = fff.readlines()
        length = len(lines)
        last = lines[-1]
        last = last.split(",")
        Emed = np.zeros(Navg)
        for k in xrange(0,Navg):
            split = lines[-1-k]
            split = split.split(",")
            Emed[k] = float(split[1])
        med = np.median(Emed)
        if med != med:
            med = 1
        try:
            ET[j] = float(elapsed[-2])/3600
            dE[j] = abs(med)
            a[j] = float(last[2])
        except:
            print 'file: '+f+' will not be included in dataset'
    if y_choice == 0:
        y = ET
    elif y_choice == 1:
        y = dE
    elif y_choice == 2:
        y = a
    if x_choice == 0:
        xvals, y = sort(xvals,y)
    plt.plot(xvals, y, marker,label=labels[i])

if y_choice == 0:
    yname = 'elapsed time (hours)'
    oname = 'ET'
elif y_choice == 1:
    yname = 'dE/E(0)'
    oname = 'dE'
elif y_choice == 2:
    yname = 'final semi-major axis of planet'
    oname = 'PL'
    plt.ylim([3,7])
#plt.yscale('log')
#plt.xscale('log')
plt.ylabel(yname, fontsize = 16)
plt.xlabel(xname, fontsize = 16)
#plt.title(title)
plt.legend(loc='upper left',prop={'size':10}, numpoints=1, markerscale=2)
plt.savefig(dirP+ext[0]+'_v_'+oname+'.png')
plt.show()

