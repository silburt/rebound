#This code samples in dt/N_planetesimal space to find convergence of simulations.

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np
#import matplotlib.pyplot as plt

#Specify what runs you want *****************************

#equal spacing
t = np.logspace(-4.25,-2,6)
N = np.logspace(2.6,4.5,6)
array = []
for tt in t:
    for NN in N:
        array.append((NN,tt))

N_runs = len(array)
Npl,timestep = zip(*array)

#plt.scatter(Npl,timestep)
#plt.xlabel('N_planetesimals')
#plt.ylabel('timestep')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim([1e2,1e5])
#plt.ylim([1e-1,1e-5])
#plt.savefig('output/sensitests.png')

runs = []
for i in xrange(0,N_runs):
    dt = timestep[i]
    Np = int(np.round(Npl[i]))
    seed = "{:.0f}".format(int(1000*random.random()))
    name = 'output/sensitests_dt%.2e_Np%d_sd%s'%(timestep[i],Np,seed)
    runs.append((dt,Np,seed,name))

os.system('make')

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+' '+str(pars[3]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=N_runs)
    pool.map(execute, runs)
    pool.close()
    pool.join()



#random generator
#np.random.seed(10)
#runs = []
#timestep = np.empty(0)
#Npl = np.empty(0)
#i=0
#while i < N_runs:
#    dt = np.random.uniform(1,4.5)
#    Np = np.random.uniform(1,5)
#    if (dt - 3)**2 + (Np - 3)**2 < 4:
#        timestep = np.append(timestep, dt-6)
#        Npl = np.append(Npl, Np)
#        i += 1