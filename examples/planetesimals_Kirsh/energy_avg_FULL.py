#This macro calculates the average energy for a set of runs. Assumes that all the files are in energy_avg/.

import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cmx
import matplotlib.colors as colors
import sys
import os
import re

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def get_cmap(N):
    '''Returns a function that maps each index in 0, 1, ... N-1 to a distinct
        RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap='hsv')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color

plot_choice = 2     #1 = plot energy, 2 = plot semi-major axis
time_sort = 0        #sort runs according to time?

ms = 0.4
alpha = 0.5
naming = ['time (years)','dE/E(0)','semi-major axis of planet (AU)']
outputn = ['time (years)','dE','a']

N_files = 0
dirP = str(sys.argv[1])
files = glob.glob(dirP+'*elapsedtime.txt')
files = sorted(files, key = natural_key)
data = []
n_it = 10e10

g_c = 0
b_c = 0
r_c = 0
if time_sort == 1:
    colors = []
time_array = ['dt12.57']
names=[]
for f in files:
    try:
        ff = open(f.split('_elapsedtime.txt')[0]+'.txt', 'r')
        lines = ff.readlines()
        length = len(lines)
        if length < n_it:   #need to find array with shortest length
            n_it = length
        data.append(lines)
        N_files += 1
        split = f.split("_")
        if time_sort == 1:
            if split[-3] == time_array[0]:
                colors.append(colorsg[g_c])
                g_c += 1
            if split[-3] == time_array[1]:
                colors.append(colorsr[r_c])
                r_c += 1
            if split[-3] == time_array[2]:
                colors.append(colorsb[b_c])
                b_c += 1
        names.append(split[-2]+', '+split[-3])
    except:
        print 'couldnt read in data file '+f

D = np.zeros(shape=(N_files,n_it))
Davg = np.zeros(n_it)
time = np.zeros(n_it)
vals_for_med = np.zeros(N_files)
for i in xrange(0,n_it):
    for j in range(0,N_files):
        split = data[j][i].split(",")
        vals_for_med[j] = float(split[plot_choice])
        D[j][i] = vals_for_med[j]
    Davg[i] = np.median(vals_for_med)
    time[i] = float(split[0])/6.283185

for i in xrange(0,N_files):
    plt.plot(time,D[i], ms=ms, color='lightgreen', alpha=alpha)
plt.plot(time,Davg, color='darkgreen', label='Average')

##############################################
#Final plotting stuff
plt.legend(loc='upper right',prop={'size':10})
plt.ylabel(naming[plot_choice], fontsize=13)
plt.xlabel('time (years)', fontsize=13)
if plot_choice == 1:
    plt.yscale('log')
    plt.xscale('log')
else:
    plt.ylim([18,26])
plt.xlim([0.5,70000])
plt.savefig(dirP+'Kirsh_avg_'+outputn[plot_choice]+'.pdf')
plt.show()