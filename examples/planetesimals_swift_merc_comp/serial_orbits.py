#Compute all the orbits.py for a given folder
import os
import sys
import time
import random
import glob

#Specify directory
dir = sys.argv[1]
files = glob.glob(dir+'*sd*.txt')
N = len(files)
i=0
while i < N:            #just want the main .txt files
    f = files[i]
    string = f.split("_")
    if string[-1]=="info.txt" or string[-1]=="elapsedtime.txt" or string[-1]=="removed.txt":
        files.remove(files[i])
        N -= 1
    else:
        i += 1

print 'computing '+str(N)+' orbit.py analyses'
for f in files:
    os.system('python orbits.py '+f)

