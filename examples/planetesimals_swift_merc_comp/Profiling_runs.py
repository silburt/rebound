#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
import numpy as np
import random

#May31 Mercury_swifter comps
tmax=62.84
Np=30000
N_runs = 3
params = []
for i in range(N_runs):
    seed = int(1000*random.random())
    name = "output/t%.0e_Np%d_sd%d"%(tmax,Np,seed)
    params.append((tmax,Np,seed,name))

N_procs = 1
os.system('make')

def execute(pars):
    mercury_dir = '../../../mercury6/input_files/Np'+str(pars[1])+'_sd'+str(pars[2])+'/'
    swifter_dir = '../../../swifter/example/input_files/Np'+str(pars[1])+'_sd'+str(pars[2])+'/'
    os.system('mkdir '+mercury_dir)
    os.system('mkdir '+swifter_dir)
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+ ' '+str(pars[3])+' '+mercury_dir+' '+swifter_dir)

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    rmv_dflt = 1
    rmv = raw_input("WARNING! Do you want to remove swifter/mercury directories? (default = 1 = yes, 0 = no) ")
    if not rmv:
        input = rmv_dflt
    else:
        input = int(rmv)
    if input == 1:
        os.system('rm -rf ../../../mercury6/input_files/*')
        os.system('rm -rf ../../../swifter/example/input_files/*')
    pool = mp.Pool(processes=N_procs)
    args=[params[i] for i in xrange(0,N_runs)]
    pool.map(execute, args)
    pool.close()
    pool.join()

#params=[(tmax,40000,13,"output/t1e5_Np40000_sd13"),(tmax,40000,14,"output/t1e5_Np40000_sd14"),(tmax,40000,50,"output/t1e5_Np40000_sd50"),(tmax,40000,60,"output/t1e5_Np40000_sd60")]
