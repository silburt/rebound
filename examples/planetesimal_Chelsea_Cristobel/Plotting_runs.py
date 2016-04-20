#Testing the dependent variable - either dt or HSR as the independent variable for various runs. Make sure the appropriate settings have been put in your problem.c file too!

import multiprocessing as mp
import os
import sys
import time
import random

#Specify what runs you want *****************************
#legend:[dt]

dt = [  6.30957344e-05,   8.57695899e-05,   1.16591440e-04,
          1.58489319e-04,   2.15443469e-04,   2.92864456e-04,
          3.98107171e-04,   5.41169527e-04,   7.35642254e-04,
          1.00000000e-03,   1.35935639e-03,   1.84784980e-03,
          2.51188643e-03,   3.41454887e-03,   4.64158883e-03,
          6.30957344e-03,   8.57695899e-03,   1.16591440e-02,
          1.58489319e-02,   2.15443469e-02,   2.92864456e-02,
          3.98107171e-02,   5.41169527e-02,   7.35642254e-02,
          1.00000000e-01]
N = len(dt)
HSR = random.sample(range(0,N),N)
HSR[:] = [2*x/float(N)+1 for x in HSR]   #HSR between 1 and 3
runs = []
for i in xrange(0,N):
    dtname = "{:.2e}".format(dt[i])
    HSRname = "{:.1f}".format(HSR[i])
    name = 'output/ChelseaChristobel_dt'+dtname+'_HSR'+HSRname
    runs.append((dt[i],HSR[i],name))

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
