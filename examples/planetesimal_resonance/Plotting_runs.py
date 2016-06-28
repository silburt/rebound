#Testing how elapsed time scales with N_planetesimals

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

#Specify what runs you want *****************************
#*****************************
N_runs = 10

random.seed()
runs = []
name = 'output/June28'
for i in xrange(0,N_runs):
    seed = "{:.0f}".format(int(1000*random.random()))
    runs.append((name,seed))

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[runs[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
