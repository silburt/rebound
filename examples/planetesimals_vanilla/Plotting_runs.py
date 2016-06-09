#Testing the dependent variable - either dt or HSR as the independent variable for various runs. Make sure the appropriate settings have been put in your problem.c file too!

import multiprocessing as mp
import os
import sys
import time
import random

#Specify what runs you want *****************************

#legend:[number of runs, min powerlaw, max powerlaw]
#params = [20,-1.5,3.5]     #HSR
params = [15,-1.5,1.5]

#*****************************
random.seed()
runs = []
lin_output = (params[2] - params[1])/float(params[0])
arg1 = params[1]
for i in xrange(0,params[0]):
    arg1 += lin_output
    seed = "{:.0f}".format(int(1000*random.random()))
    name = 'output/Vanilla_Jup2edgepowerlaw'+str(arg1)+'_Np20000_HSR6_Mpl10mp'
    runs.append((str(arg1),seed,name))

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
