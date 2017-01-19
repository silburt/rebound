#Testing the dependent variable - either dt or HSR as the independent variable for various runs. Make sure the appropriate settings have been put in your problem.c file too!

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

#Specify what runs you want *****************************

N_runs = 200
HSF = 3
theta = np.linspace(0,N_runs)
f = random.sample(np.linspace(0,1,N_runs), N_runs)

runs = []
dir = 'output/'
for i in xrange(0,N_runs):
    name = dir+'HSF3_theta%.2f_f%.2f'%(theta[i],f[i])
    runs.append([HSF,name,theta[i],f[i]])

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+ ' '+str(pars[3]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=25)
    args=[runs[i] for i in xrange(0,N_runs)]
    pool.map(execute, args)
    pool.close()
    pool.join()
