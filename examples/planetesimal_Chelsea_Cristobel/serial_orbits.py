#Compute all the orbits.py for a given folder

import multiprocessing as mp
import os
import sys
import time
import random
import glob

ext = ''

#Specify directory
dir = sys.argv[1]
files = glob.glob(dir+'*'+ext+'.txt')
N = len(files)
i=0
while i < N:    #just want the main .txt files
    f = files[i]
    string = f.split("_")
    if string[-1]=="removed.txt" or string[-1]=="elapsedtime.txt":
        files.remove(files[i])
        N -= 1
    else:
        i += 1

def execute(filename):
    os.system('python orbits.py '+filename+' 1')
    print 'finished process'

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=N)
    print 'computing '+str(N)+' orbit.py analyses'
    pool.map(execute, files)
    pool.close()
    pool.join()