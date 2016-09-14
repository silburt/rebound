#Testing how elapsed time scales with N_planetesimals

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

#Specify what runs you want *****************************
N_K_runs = 4    #number of different K values, as in Lee & Peale (2002)
N_a_runs = 6    #number of different super-Earth positions

random.seed()
runs = []
K_vals = np.round(np.logspace(0,2,N_K_runs))
K = np.empty(0)
for i in range(N_a_runs):
    K = np.concatenate((K,K_vals))

a_vals = np.linspace(3,4.5,N_a_runs)
a = np.empty(0)
for i in range(N_a_runs):
    a = np.concatenate((a,np.ones(N_K_runs)*a_vals[i]))

name = 'output/SuperEarth'
for i in xrange(0,N_K_runs*N_a_runs):
    name = 'output/SuperEarth_asuper%.1f_K%d'%(a[i],K[i])
    runs.append((a[i],K[i],name))

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    pool.map(execute, runs)
    pool.close()
    pool.join()
