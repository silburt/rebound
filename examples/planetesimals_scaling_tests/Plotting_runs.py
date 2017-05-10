#Testing the dependent variable - either dt or HSR as the independent variable for various runs. Make sure the appropriate settings have been put in your problem.c file too!

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

#Specify what runs you want *****************************

N_runs = 50
HSF = np.logspace(-2,1.7,N_runs)
theta = np.linspace(0,1,N_runs)
ff = random.sample(np.linspace(0,1,N_runs), N_runs)

runs = []
dir = 'output/'
for i in xrange(0,N_runs):
    name = dir+'HSF%.1e_theta%.2f_f%.2f'%(HSF[i],theta[i],ff[i])
    runs.append([HSF[i],name,theta[i],ff[i]])

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+ ' '+str(pars[3]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=5)
    args=[runs[i] for i in xrange(0,N_runs)]
    pool.map(execute, args)
    pool.close()
    pool.join()
