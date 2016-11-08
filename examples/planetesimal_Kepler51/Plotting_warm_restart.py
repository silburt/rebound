#Testing how elapsed time scales with N_planetesimals

import multiprocessing as mp
import os
import sys
import glob
import random
import numpy as np

#Specify what runs you want *****************************
runs = []

dir = 'saved_output/round1_coldvwarm/'
files = np.array(glob.glob(dir+'*t=665319.bin'))

#grab 6 random runs
length = 6
runs = files[np.random.randint(0,len(files),length)]

os.system('make')

def execute(pars):
    ext = pars.split(".bin")[0]
    os.system('./rebound '+ext)

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    pool.map(execute, runs)
    pool.close()
    pool.join()

