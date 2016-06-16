#Testing how elapsed time scales with N_planetesimals

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

#Specify what runs you want *****************************
#*****************************
random.seed()
runs = []
Np = np.logspace(1,4,50,dtype=int)
for i in xrange(0,len(Np)):
    seed = "{:.0f}".format(int(1000*random.random()))
    name = 'output/Np'+str(Np[i])+'_HSR3_dt0.05'
    runs.append((Np[i],seed,name))

os.system('make')

length = len(runs)

def execute(pars):
    print str(pars[0]), str(pars[1]), str(pars[2])
    #os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    args=[runs[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
