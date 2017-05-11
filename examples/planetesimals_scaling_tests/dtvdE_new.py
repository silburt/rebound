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
    coeff = (dt*dt / 12)*(Me*mp*Ms/(rhHSR**2))
    #theory
    term1 = 1./(a*rhHSR)
    term2 = Me/(Ms*rhHSR**2)
    term3 = 1./(2*a*a)
    theory = (term1 + term2 + term3)*coeff*scale_factor
    #plot
    if choice == 0:
        x = dt
    else:
        x = HSR
    plt.plot(x, theory,label='theory',linewidth=3)

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

#plot dt vs.
dir = sys.argv[1]
y_choice = int(sys.argv[2])      #0 = plot elapsed time, 1 = energy
x_choice = 1                     #0 = dt, 1 = HSR, 2 = Np
Navg = 3                        #number of points at end of .txt file to average energy over

files = glob.glob('%s*elapsedtime.txt'%dir)

a = np.array(1)             #AU
Ms = np.array(1)            #Solar mass
mp = np.array(1e-8)         #Solar mass
Me = np.array(5e-5)         #Solar mass
rh = a*(Me/(3*Ms))**(1./3.) #AU

dE, HSF, dt, ET = [], [], [], []
for f in files:
    fd=f.split('_elapsedtime.txt')[0] + '.txt'
    try:
        data = np.genfromtxt(fd,delimiter=',')
        dE.append(np.mean(data[-Navg-1:-1,1]))
        HSF.append(data[0,3])
        dt.append(data[0,4])
        ET.append(float(open(f,'r').readlines()[1].split()[-2]))
    except:
        print 'error processing %s'%fd

if y_choice == 0:
    yname = 'elapsed time (seconds)'
    yvar = ET
elif y_choice == 1:
    name = 'relative energy error'
    yvar = dE

if x_choice == 0:
    xname = 'timestep (yr/2pi)'
    xvar = dt
elif x_choice == 1:
    xname = 'HSF (Hill Radii)'
    xvar = HSF

plt.plot(xvar, yvar, '.')
plt.yscale('log')
plt.xscale('log')
plt.ylabel(name, fontsize = 13)
plt.xlabel(xname, fontsize = 13)
#plt.title(title)
#plt.ylim([np.min(y),2])
plt.legend(loc='upper left',prop={'size':10}, numpoints=1, markerscale=2)
plt.savefig('%sdtvdE_new.pdf'%dir)
plt.show()


'''
#old
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
    '''
'''
    #older
def theory(a,Ms,mp,Me,dt,rh,HSR,choice,scale_factor):
    print scale_factor
    rhHSR = rh*HSR
    p = mp/np.sqrt(a)
    a2 = a*a
    ar = (a*a - rhHSR*rhHSR)**0.5
    M3 = Me*mp*Ms
    tau2 = dt*dt / 12.
    #theory
    original = Me*p*p/(mp*rhHSR**3) + 0.5*Me*Me*mp/rhHSR**4 - 0.5*Me*mp/(rhHSR*ar**3)
    term1 = -(0.5*p*p/mp)*(1./(ar**3) + 3*rhHSR*rhHSR/ar**5)
    term2 = mp*rhHSR*rhHSR/(2*ar**6)
    term3 = 0.5*Me*mp/rhHSR
    theory = (original + term3)*tau2*scale_factor
    #theory = (term1 + term2+ term3)*tau2*scale_factor
    #plot
    if choice == 0:
    x = dt
    else:
    x = HSR
    plt.plot(x, theory,'+-',label='R3B theory')
    '''
