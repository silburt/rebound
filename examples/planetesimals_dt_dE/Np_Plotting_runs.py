#Testing how energy error scales with Np! Use in HSR mode.

import multiprocessing as mp
import os
import sys
import time
import random

#Specify what runs you want *****************************

#legend:[HSR,N_planetesiamls_max, arg1 Type, # phase changes]
HSR = 3
Np_ini = 1     #HSR
Np_fini = 1000
N_runs = 20
log_output = (Np_fini/Np_ini)**(1./N_runs)
#*****************************
seed = "12"
runs = []
Np = Np_ini
for i in xrange(0,N_runs):
    #theta = "{:.2f}".format(random.random()) #range from 0-1, only matters for the single planetesimal tests
    theta = "0"
    name = 'output/Np'+str(Np)+'_HSR'+str(HSR)+'_th'+theta
    runs.append((HSR,Np,seed,name,theta))
    Np = int(Np*log_output + 1)

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+ ' '+str(pars[3])+' '+str(pars[4]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[runs[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
