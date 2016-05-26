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

y_choice = int(sys.argv[2])      #1 = energy, 2=energy / sqrt(N_CE)

dirP = sys.argv[1]
Navg = 3                       #number of points at end of .txt file to average energy over
xvals = []

a = np.array(1)             #AU
Ms = np.array(1)            #Solar mass
mp = np.array(2.5e-7)         #Solar mass
Me = np.array(5e-5)         #Solar mass
rh = a*(Me/(3*Ms))**(1./3.) #AU
HSR = 6
dt = 0.015/twopi

files = glob.glob(dirP+'*_elapsedtime.txt')
lenf = len(files)
ET = np.zeros(lenf)
dE = np.zeros(lenf)
CE = np.zeros(lenf)
ratio = np.zeros(lenf)  #ratio of energy to N_CE
scaling_factor = np.zeros(0)
for j,f in enumerate(files):
    fff = open(f.split("_elapsedtime.txt")[0]+".txt", 'r')
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
    med = np.median(Emed)
    if med != med:
        med = 1
    try:
        dE[j] = med
        CE[j] = float(last[5])
        if CE[j] == 0:
            CE[j] = 1
        ratio[j] = dE[j]/np.sqrt(CE[j])
            #ratio[j] = dE[j]/CE[j]
    except:
        print 'file: '+f+' will not be included in dataset'

if y_choice == 0:
    y = dE
    yname = 'Fractional Energy Error'
elif y_choice == 1:
    y = ratio
    yname = 'dE/E/N$_{switch}$'
plt.plot(CE, y, '.')

#theory
rhHSR = rh*HSR
M3 = Me*mp*Ms
tau2 = dt*dt / 12.
term1 = M3/(a*rhHSR**3)
term2 = Me*Me*mp/(rhHSR**4)
term3 = -M3/(2*rhHSR*((a*a - rhHSR**2)**1.5))    #minor term, neg^x returns invalue value
termp = 3*M3/((a*rhHSR**3))
theory = (term1 + term2 + term3)*tau2*np.mean(scaling_factor)
plt.plot(CE, theory*np.sqrt(CE),'+-',label='$E^{WH}_{floor} * \sqrt{CE}$')
plt.plot(CE, theory*CE,'<-',markeredgecolor='none',label='$E^{WH}_{floor} * CE$')

plt.yscale('log')
plt.xscale('log')
plt.ylabel(yname, fontsize = 13)
plt.xlabel('Number of Close Encounters', fontsize = 13)
plt.xlim([5,1e4])
plt.legend(loc='upper left',prop={'size':13}, numpoints=1, markerscale=2)
#plt.title('HSR='+str(HSR)+', dt='+str(round(dt*twopi,4))+', a$_p$='+str(a)+', N$_{pl}$=200')
plt.savefig(dirP+'CEvdE.pdf')
plt.show()

