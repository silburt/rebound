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

def theory(a,Ms,mp,Me,dt,rh,HSR,choice,scale_factor):
    print scale_factor
    rhHSR = rh*HSR
    a2 = a*a
    M3 = Me*mp*Ms
    tau2 = dt*dt / 12.
    #theory
    term1 = M3/(a*rhHSR**3)
    term2 = Me*Me*mp/(rhHSR**4)
    term3 = -M3/(2*rhHSR*((a2 - rhHSR**2)**1.5))    #minor term, neg^x returns invalue value
    termp = 3*M3/((a*rhHSR**3))
    theory = (term1 + term2 + term3)*tau2*scale_factor
    #theory = (termp)*tau2*scale_factor
    #plot
    if choice == 0:
        x = dt
    else:
        x = HSR
    plt.plot(x, theory,'+-',label='R3B theory')

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

#plot dt vs.
y_choice = int(sys.argv[1])      #0 = plot elapsed time, 1 = energy, 2=N_CE
x_choice = 1                     #0 = dt, 1 = HSR, 2 = Np

dirP = str('dtvdE_files/')
Navg = 3                       #number of points at end of .txt file to average energy over
xvals = []

if x_choice == 0:   #dt
    ext = 'dt'
    xname = 'timestep (years)'
    yvals = ['Np0']
    path = '_th*_elapsedtime.txt'
    labels = ['Simulation']
    marker = '.'
    title = 'dt: Integrating 2p, 2pl system for 1000 orbits'
    HSR = 6
if x_choice == 1:   #HSR
    ext = 'HSR'
    xname = 'HSB (hill radii)'
    #yvals = ['Np500','Np5000']
    #path = '_th*_elapsedtime.txt'
    yvals = ['Np1000_sd*']
    #yvals = ['Np0_th*']
    path = '_elapsedtime.txt'
    labels = ['Simulation']
    marker = '.'
    title = 'HSR: Integrating 2p, 2pl system for 1M orbits'
    dt = np.array(1e-3)/twopi  #must be in yr/2pi format
if x_choice == 2:   #Np
    ext = 'Np'
    xname = 'Number of Planetesimals'
    yvals = ['']
    path = '_HSR6.00_th*_elapsedtime.txt'
    labels = yvals
    marker = '.'
    title = 'dt: Integrating 2p, 2pl system for 1000 orbits'
    HSR = 6
    dt = np.array(1e-3)/twopi #must be in yr/2pi format

a = np.array(1)             #AU
Ms = np.array(1)            #Solar mass
mp = np.array(1e-8)         #Solar mass
Me = np.array(5e-5)         #Solar mass
rh = a*(Me/(3*Ms))**(1./3.) #AU

leny = len(yvals)
for i in xrange(0,leny):
    files = glob.glob(dirP+'*'+yvals[i]+path)
    files = sorted(files, key = natural_key)    #assuming files are sorted in correct order now
    lenf = len(files)
    ET = np.zeros(lenf)
    dE = np.zeros(lenf)
    CE = np.zeros(lenf)
    ratio = np.zeros(lenf)  #ratio of energy to N_CE
    xvals = np.zeros(lenf)
    scaling_factor = np.zeros(0)
    for j in xrange(0,lenf):
        f = files[j]
        header = f.split("_")
        xx = header[-4]
        xvals[j] = float(re.sub('^files/'+ext, '', xx))
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
        if float(last[6]) != 0:
            scaling_factor = np.append(scaling_factor,float(last[1])/float(last[6]))
        Emed = np.zeros(Navg)
        for k in xrange(0,Navg):
            split = lines[-1-k]
            split = split.split(",")
            Emed[k] = float(split[1])
            #Emed[k] = float(split[6])
        med = np.median(Emed)
        if med != med:
            med = 1
        try:
            ET[j] = float(elapsed[-2])/3600
            dE[j] = med
            CE[j] = float(last[5])
            if x_choice != 2:
                if CE[j] == 0:
                    CE[j] = 1
                #ratio[j] = dE[j]/np.sqrt(CE[j])
                ratio[j] = dE[j]/CE[j]
        except:
            print 'file: '+f+' will not be included in dataset'
    if y_choice == 0:
        y = ET
    elif y_choice == 1:
        y = dE
    elif y_choice == 2:
        y = CE
    else:
        y = ratio
    if x_choice == 0:
        xvals, y = sort(xvals,y)
    plt.plot(xvals, y, marker,label=labels[i])
if x_choice == 1 and y_choice == 1 or y_choice == 3:
    theory(a,Ms,mp,Me,dt,rh,xvals,x_choice,np.mean(scaling_factor))
elif x_choice == 0 and y_choice == 1 or y_choice == 3:
    #xvals /= twopi
    theory(a,Ms,mp,Me,xvals,rh,HSR,x_choice,np.mean(scaling_factor))

if y_choice == 0:
    name = 'elapsed time (seconds)'
    oname = 'ET'
elif y_choice == 1:
    name = 'Fractional Energy Error'
    oname = 'dE'
elif y_choice == 2:
    name = '# Close Encouters'
    oname = 'CE'
else:
    name = 'dE/E/N$_{switch}$'
    oname = 'ratio'
plt.yscale('log')
plt.xscale('log')
plt.ylabel(name, fontsize = 13)
plt.xlabel(xname, fontsize = 13)
#plt.title(title)
#plt.ylim([np.min(y),2])
plt.legend(loc='upper left',prop={'size':10}, numpoints=1, markerscale=2)
plt.savefig('dtvdE_files/'+ext+'_v_'+oname+'.pdf')
plt.show()

