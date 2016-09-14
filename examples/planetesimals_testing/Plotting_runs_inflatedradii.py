#Testing the dependent variable - either dt or HSR as the independent variable for various runs. Make sure the appropriate settings have been put in your problem.c file too!

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

#Specify what runs you want *****************************
n_runs = 2
N_planetesimals = 1000
tmax = 1e6
radius = np.logspace(0,1.8,n_runs)
radius = np.insert(radius,1,[1,1])
n_runs += 2
runs = []

for i in range(n_runs):
    seed = int(1000*np.random.random())
    name = 'output/HERMES_Np'+str(N_planetesimals)+'_rpfac%.1f'%radius[i]+'_sd'+str(seed)
    runs.append((tmax,N_planetesimals,radius[i],seed,name))

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+' '+str(pars[3])+' '+str(pars[4]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=n_runs)
    pool.map(execute, runs)
    pool.close()
    pool.join()