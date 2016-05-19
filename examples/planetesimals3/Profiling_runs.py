#Testing how energy and time scales with N_particles.

import multiprocessing as mp
import os
import sys
import time
#import numpy as np

#params=[(10000000,50,12,"output/Np50_sd12.txt"),(10000000,50,22,"output/Np50_sd22.txt"),(10000000,50,32,"output/Np50_sd32.txt"),(10000000,50,51,"output/Np50_sd51.txt"),(10000000,50,61,"output/Np50_sd61.txt"),(10000000,50,80,"output/Np50_sd80.txt"),(5000000,100,10,"output/Np100_sd10.txt"),(5000000,100,18,"output/Np100_sd18.txt"),(5000000,100,27,"output/Np100_sd27.txt"),(5000000,100,41,"output/Np100_sd41.txt"),(5000000,100,73,"output/Np100_sd73.txt"),(5000000,100,83,"output/Np100_sd83.txt"),(1000000,500,13,"output/Np500_sd13.txt"),(1000000,500,16,"output/Np500_sd16.txt"),(1000000,500,33,"output/Np500_sd33.txt"),(1000000,500,39,"output/Np500_sd39.txt"),(1000000,500,79,"output/Np500_sd79.txt"),(1000000,500,90,"output/Np500_sd90.txt")]

#Mercury/Hybarid/Swifter comp tests
#params=[(50000000,50,12,"output/Np50_sd12"),(50000000,50,22,"output/Np50_sd22"),(50000000,50,32,"output/Np50_sd32"),(50000000,50,51,"output/Np50_sd51"),(50000000,50,61,"output/Np50_sd61"),(50000000,50,80,"output/Np50_sd80")]

#May19 Mercury_swifter comps
tmax=5e7
params=[(tmax,50,12,"output/t5e7_Np50_sd12"),(tmax,50,22,"output/t5e7_Np50_sd22"),(tmax,50,32,"output/t5e7_Np50_sd32"),(tmax,50,51,"output/t5e7_Np50_sd51"),(tmax,50,61,"output/t5e7_Np50_sd61"),(tmax,50,80,"output/t5e7_Np50_sd80")]

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
