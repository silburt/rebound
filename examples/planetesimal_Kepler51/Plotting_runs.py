#Testing how elapsed time scales with N_planetesimals

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

#Specify what runs you want *****************************
N_runs = 6

random.seed()
runs = []
Np = "{:.0f}".format(2000)
for i in xrange(0,N_runs):
    seed = "{:.0f}".format(int(1000*random.random()))
    name = 'output/cold_Np%s_sd%s'%(Np,seed)
    runs.append((0,Np,seed,name))

for i in xrange(0,N_runs):
    seed = "{:.0f}".format(int(1000*random.random()))
    name = 'output/warm_Np%s_sd%s'%(Np,seed)
    runs.append((1,Np,seed,name))

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+' '+str(pars[3]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    pool.map(execute, runs)
    pool.close()
    pool.join()
