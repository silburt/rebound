#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
import numpy as np
import random

#May31 Mercury_swifter comps
#tmax=5e7
#params=[(tmax,100,13,"output/t5e7_Np100_sd13"),(tmax,100,23,"output/t5e7_Np100_sd23"),(tmax,100,33,"output/t5e7_Np100_sd33"),(tmax,100,52,"output/t5e7_Np100_sd52"),(tmax,100,62,"output/t5e7_Np100_sd62"),(tmax,100,85,"output/t5e7_Np100_sd85")]

params = []
np.random.seed()
tmax = "{:.0e}".format(10000)
Np = np.logspace(1,4,50,dtype=int)
for i in xrange(0,len(Np)):
    seed = "{:.0f}".format(int(1000*random.random()))
    name = 'output/t'+tmax+'_Np'+str(Np[i])+'_HSR3_dt0.01'
    params.append((tmax,Np[i],seed,name))

length = len(params)

os.system('make')

def execute(pars):
    mercury_dir = '../../../mercury6/input_files/Np'+str(pars[1])+'_sd'+str(pars[2])+'/'
    swifter_dir = '../../../swifter/example/input_files/Np'+str(pars[1])+'_sd'+str(pars[2])+'/'
    os.system('mkdir '+mercury_dir)
    os.system('mkdir '+swifter_dir)
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+ ' '+str(pars[3])+' '+mercury_dir+' '+swifter_dir)

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    rmv_dflt = 1
    rmv = raw_input("WARNING! Do you want to remove swifter/mercury directories? (default = 1 = yes, 0 = no) ")
    if not rmv:
        input = rmv_dflt
    else:
        input = int(rmv)
    if input == 1:
        os.system('rm -rf ../../../mercury6/input_files/*')
        os.system('rm -rf ../../../swifter/example/input_files/*')
    pool = mp.Pool(processes=length)
    args=[params[i] for i in xrange(0,length)]
    pool.map(execute, args)
    pool.close()
    pool.join()
