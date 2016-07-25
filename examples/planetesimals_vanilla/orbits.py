import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re
 
def cdf(array):
    return np.arange(1,len(array)+1)/float(len(array))

def removed_particles(file_name):
    split1 = file_name.split(".txt")
    ff = open(split1[0]+"_removed.txt", 'r')
    lines = ff.readlines()
    eject = np.zeros(0)
    collide = np.zeros(0)
    time = np.zeros(0)
    Ntot=0
    Ncoll=0
    Nej=0
    for line in lines:
        split = line.split(',')
        val = float(split[1])/6.283185
        if split[0] == "Collision":
            collide = np.append(collide, val)
            Ncoll += 1
        elif split[0] == "Ejection":
            eject = np.append(eject, val)
            Nej += 1
        time = np.append(time, val)
        Ntot += 1
    return time, cdf(time)*Ntot, collide, cdf(collide)*Ncoll, eject, cdf(eject)*Nej

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])
if len(sys.argv) > 2:
    suppress_plot = int(sys.argv[2])
else:
    suppress_plot = 0

diagnostics = 1
plot_dr_or_com = 1  #0 = dr, 1 = com

fos = open(file_name, 'r')
data = np.loadtxt(fos, delimiter=',')
timex, timey, collidex, collidey, ejectx, ejecty = removed_particles(file_name)

signed_energy = 1
reverse_xaxis = 0

ms=3
if diagnostics == 1:
    time,dE,ainner,Nglob,Nmini,angM,com,rinner,router,aouter = data[:,0]/6.283,data[:,1],data[:,2],data[:,3],data[:,4],data[:,6],data[:,8],data[:,9],data[:,10],data[:,11],
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,10))
    if signed_energy == 1:
        axes[0].plot(time[dE>0],abs(dE[dE>0]),'o', ms=ms, markeredgecolor='none',color='blue',label='positive energy')
        axes[0].plot(time[dE<0],abs(dE[dE<0]),'o', ms=ms, markeredgecolor='none',color='green',label='negative energy')
        axes[0].plot(time[angM>0],abs(angM[angM>0]),'o', ms=ms, markeredgecolor='none',color='red',label='positive ang. mom.')
        axes[0].plot(time[angM<0],abs(angM[angM<0]),'o', ms=ms, markeredgecolor='none',color='orange',label='negative ang. mom.')
        axes[0].set_xscale('log')
    else:
        axes[0].plot(time,abs(dE), 'o', ms=ms, markeredgecolor='none',color='blue',label='dE/E(0)')
        axes[0].plot(time,abs(angM), 'o', ms=ms, markeredgecolor='none',color='green',label='dL/L(0)')
        axes[0].plot(time,0.5e-12*time, color='red', label='t')
        axes[0].plot(time,1e-12*time**0.5, color='black', label='t^0.5')
        axes[0].set_xscale('log')
    axes[0].set_xlim([0.1,max(time)])
    axes[0].set_yscale('log')
    axes[0].legend(loc='upper left',prop={'size':10})
    axes[1].plot(timex,timey, color='black', label='total removed particles, N='+str(len(timex)))
    axes[1].plot(collidex,collidey, color='red', label='collided particles, N='+str(len(collidex)))
    axes[1].plot(ejectx,ejecty, color='blue', label='ejected particles, N='+str(len(ejectx)))
    axes[1].legend(loc='upper left',prop={'size':10})
    axes[1].set_ylabel('Removed particles')
    axes[1].set_xscale('log')
    axes[1].set_yscale('log')
    axes[1].set_xlim([0.1,max(time)])
    axes[2].plot(time,aouter, 'o', ms=ms, markeredgecolor='none')
    axes[2].set_ylabel('semi-major axis')
#com
    #axes[2].plot(time,data[:,11], 'o', ms=ms, markeredgecolor='none',label='')
    #axes[2].plot(time,1e-11*time**1.5, color='black',label='t^1.5')
    #axes[2].plot(time,1e-11*time**2, color='red',label='t^2')
    #axes[2].set_xlabel('simulation time (yrs)')
    #axes[2].set_ylabel('dcom')
    #axes[2].set_yscale('log')
    #axes[2].legend(loc='lower left',prop={'size':10})
else:
    plt.plot(data[:,0],data[:,1], 'o', ms=ms, markeredgecolor='none')
    #plt.plot(data[:,0],3e-10*data[:,0]**(0.5),color='black',label='t^1/2 growth')
    plt.ylabel('Energy')
    plt.xlabel('time (years)')
    #plt.legend(loc='upper left',prop={'size':10})
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([0.1,data[-1,0]])

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
if suppress_plot != 1:
    plt.show()