#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
import time
import numpy as np

n_runs = 6
#HSR = 4.
#dt = 4.
#rr = 7.88215e-5
#rr = 5e-4
#ones = np.ones(6)
#params = zip(ones*HSR, ones*dt, np.random.randint(0,1000,n_runs), ones*rr)
params = np.random.randint(0,1000,n_runs)

length = len(params)

os.system('make')

def execute(pars):
    #os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+' '+str(pars[3]))
    os.system('./rebound '+str(pars[0]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=n_runs)
    pool.map(execute, params)
    pool.close()
    pool.join()
