#Testing the dependent variable - either dt or HSR as the independent variable for various runs. Make sure the appropriate settings have been put in your problem.c file too!

import multiprocessing as mp
import os
import sys
import time
import random

dt = [  6.30957344e-05,   1.00000000e-04,   1.58489319e-04,
      2.51188643e-04,   3.98107171e-04,   6.30957344e-04,
      1.00000000e-03,   1.58489319e-03,   2.51188643e-03,
      3.98107171e-03,   6.30957344e-03,   1.00000000e-02]
HSR = [ 2.71608978,  2.37862452,  2.60494698,  1.1539763 ,  2.68925867,
       1.22855337,  2.56077752,  2.43665146,  2.44957389,  2.49497113,
       1.53984808,  2.14060528]
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
