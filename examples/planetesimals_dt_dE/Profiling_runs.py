#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
import time
import numpy as np
#import numpy as np

n_runs = 200
time = np.logspace(1,4,n_runs)
seed = np.random.randint(0,1000,n_runs)
names = []
for i in xrange(0,n_runs):
    names.append("output/tmax"+str(round(time[i],2))+"_Np200_HSR6_dt0.015")

params = zip(time,seed,names)
length = len(params)

os.system('make')

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    pool.map(execute, params)
    pool.close()
    pool.join()
