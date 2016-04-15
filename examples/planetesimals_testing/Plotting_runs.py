#Testing the dependent variable - either dt or HSR as the independent variable for various runs. Make sure the appropriate settings have been put in your problem.c file too!

import multiprocessing as mp
import os
import sys
import time
import random

#Specify what runs you want *****************************
'''
#legend:[number of runs, N_planetesiamls]
params = [6,500]

#*****************************
random.seed()
runs = []
Np = str(params[1])
for i in xrange(0,params[0]):
    seed = "{:.0f}".format(int(1000*random.random()))
    name = 'output/HYBARID_Np'+Np+'_sd'+seed
    runs.append((Np,seed,name))
'''

#legend:[planet mass fac]
params = range(1,10,1)
seed = random.sample(range(0,1000),len(params))
runs = []
for i in xrange(0,len(params)):
    name = 'output/HYBARID_Np500_sd'+str(seed[i])+'_mfac'+str(params[i])
    runs.append((params[i],seed[i],name))

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
