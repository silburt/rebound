#Compute all the orbits.py for a given folder

import multiprocessing as mp
import os
import sys
import time
import random
import glob

#ext = '_Mpl10mp'    #what goes right before the .txt
ext = 'sd*'

#Specify directory
dir = sys.argv[1]
Nplanets = sys.argv[2]
files = glob.glob(dir+'*'+ext+'.txt')
N = len(files)
i=0
while i < N:            #just want the main .txt files
    f = files[i]
    string = f.split("_")
    if string[-1]=="info.txt" or string[-1]=="elapsedtime.txt" or string[-2]=="eiasnapshot":
        files.remove(files[i])
        N -= 1
    else:
        i += 1

vals = []
for i in range(N):
    vals.append((files[i],Nplanets))

def execute(pars):
    #os.system('python orbits.py '+pars[0]+' '+pars[1])
    os.system('python check_resonance.py '+pars[0]+' '+pars[1])
    print 'finished process'

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=N)
    print 'computing '+str(N)+' orbit.py analyses'
    pool.map(execute, vals)
    pool.close()
    pool.join()
