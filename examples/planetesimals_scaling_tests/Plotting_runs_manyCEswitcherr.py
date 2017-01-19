#Testing how elapsed time scales with N_planetesimals

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

#Specify what runs you want *****************************
#*****************************
random.seed()
runs = []
N_runs = 4
tmax = 1e6
for i in xrange(0,N_runs):
    seed = int(1000*random.random())
    name = 'output/t%.0e_AutoHSF3_dt0.02_sd%d'%(tmax, seed)
    runs.append((tmax,seed,name))

os.system('make')

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=N_runs)
    args=[runs[i] for i in xrange(0,N_runs)]
    pool.map(execute, args)
    pool.close()
    pool.join()
