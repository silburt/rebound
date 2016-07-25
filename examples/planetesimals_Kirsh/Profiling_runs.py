#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
import time
import numpy as np

#params=[(4,10,0.25),(4,10,0.5),(4,10,1),(4,10,2),(4,10,4),(4,10,8),(4,10,16),(4,20,2),(4,30,2),(4,40,2),(4,50,2),(4,60,2),(12.76,10,0.25),(12.76,10,0.5),(12.76,10,1),(12.76,10,2),(12.76,10,4),(12.76,10,8),(12.76,10,16),(12.76,20,2),(12.76,30,2),(12.76,40,2),(12.76,50,2),(12.76,60,2)]

n_runs = 6
HSR = 4.
dt = 4.
rr = 7.88215e-5
#rr = 5e-4

ones = np.ones(6)
params = zip(ones*HSR, ones*dt, np.random.randint(0,1000,n_runs), ones*rr)

length = len(params)

os.system('make')

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+' '+str(pars[3]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=n_runs)
    pool.map(execute, params)
    pool.close()
    pool.join()