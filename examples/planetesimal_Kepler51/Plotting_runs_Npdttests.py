#This code samples in dt/N_planetesimal space to find convergence of simulations.

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np
import matplotlib.pyplot as plt

#Specify what runs you want *****************************
N_runs = 40
np.random.seed(10)
runs = []
timestep = np.empty(0)
Npl = np.empty(0)
i=0
while i < N_runs:
    dt = np.random.uniform(1,4.5)
    Np = np.random.uniform(1,5)
    if (dt - 3)**2 + (Np - 3)**2 < 4:
        timestep = np.append(timestep, dt-6)
        Npl = np.append(Npl, Np)
        i += 1

plt.scatter(Npl,timestep)
plt.show()

for i in xrange(0,N_runs):
    dt = "{:.1e}".format(10**(timestep[i]))
    Np = "{:.0e}".format(10**(Npl[i]))
    seed = "{:.0f}".format(int(1000*random.random()))
    name = 'output/cold_dt%s_Np%s_sd%s'%(dt,Np,seed)
    runs.append((dt,Np,seed,name))

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+' '+str(pars[3])+' '+str(pars[4]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    pool.map(execute, runs)
    pool.close()
    pool.join()
