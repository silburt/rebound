#Testing the dependent variable - either dt or HSR as the independent variable for various runs. Make sure the appropriate settings have been put in your problem.c file too!

import multiprocessing as mp
import os
import sys
import time
import random

dt = [1e-3,1e-4,1e-3,1e-4,1e-3,1e-4,1e-3,1e-4,1e-3,1e-4,1e-3,1e-4]
HSR = [1,1,1.05,1.05,1.1,1.1,1.15,1.15,1.2,1.2,1.25,1.25]
runs = []
for i in xrange(0,len(dt)):
    dtstr = "{:.4f}".format(dt[i])
    HSRstr = "{:.2f}".format(HSR[i])
    name = 'output/Chelsea_HSR'+HSRstr+'_dt'+dtstr
    runs.append((HSRstr,dtstr,name))

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[runs[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
