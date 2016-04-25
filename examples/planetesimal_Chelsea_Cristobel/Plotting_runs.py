#Testing the dependent variable - either dt or HSR as the independent variable for various runs. Make sure the appropriate settings have been put in your problem.c file too!

import multiprocessing as mp
import os
import sys
import time
import random

dt = [  1.00000000e-06,   2.15443469e-06,   4.64158883e-06,
      1.00000000e-05,   2.15443469e-05,   4.64158883e-05,
      1.00000000e-04,   2.15443469e-04,   4.64158883e-04,
      1.00000000e-03]
HSR = [ 3.99436153,  3.29998598,  3.00573604,  3.47210691,  3.90949609,
       3.95028692,  3.2805955 ,  2.17323678,  2.64272983,  3.48236786]
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
