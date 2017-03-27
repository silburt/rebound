#Testing the dependent variable - either dt or HSR as the independent variable for various runs. Make sure the appropriate settings have been put in your problem.c file too!

import multiprocessing as mp
import os
import sys
import time
import random

#Specify what runs you want *****************************

#legend:[number of runs, tmax,  N_planetesiamls]
params = [4,5e7,50]

#*****************************
random.seed()
runs = []
tmax = "{:.0e}".format(params[1])
Np = str(params[2])
for i in xrange(0,params[0]):
    seed = "{:.0f}".format(int(1000*random.random()))
    name = 'output/t'+tmax+'_Np'+Np+'_sd'+seed
    runs.append((tmax,Np,seed,name))

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+ ' '+str(pars[3]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[runs[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
