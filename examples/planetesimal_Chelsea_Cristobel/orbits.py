import sys
import matplotlib.pyplot as plt
import numpy as np
import re

plot = int(sys.argv[1])

if plot == 0:
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,10), sharex=True)
    file_name = "output/test.txt"
    fos = open(file_name, 'r')
    time, dE, N_mini, N, ET, HSF = np.loadtxt(fos, delimiter=',', unpack=True)
    
    axes[0].plot(time,dE, '.')
    axes[0].plot(time,0.5e-10*time, color='red', label='t')
    axes[0].plot(time,2e-9*time**0.5, color='black', label='t^0.5')
    axes[0].set_xlim([0.1,max(time)])
    axes[0].set_yscale('log')
    axes[1].plot(time,HSF, '.')
    axes[1].set_ylabel('HSF')
    axes[2].plot(time,N,'.',label='N')
    axes[2].plot(time,N_mini,'.',label='N_mini')
    axes[2].legend(loc='upper left')
    axes[0].set_xscale('log')

else:
    fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(10,10), sharex=True)
    file_name = "output/test_HSF.txt"
    fos = open(file_name, 'r')
    time, min_dt_ratio, HSF, vmax, hillmax, c1,c2,c3,c4 = np.loadtxt(fos, delimiter=',', unpack=True)

    axes[0].plot(time,HSF, '.')
    axes[0].set_ylabel('HSF')
    axes[1].plot(time,min_dt_ratio, label='min_dt_enc/r->dt')
    axes[1].plot([0,max(time)],[1,1],'r--', label='threshold for HSF increase')
    axes[1].legend(loc='upper left',fontsize=10)
    axes[1].set_yscale('log')
    axes[1].set_ylim([0.5*min(min_dt_ratio),2*max(min_dt_ratio)])
    axes[2].plot(time,vmax, '.', label='vmax')
    axes[2].plot(time,hillmax, '.', label='hillmax')
    axes[2].legend(loc='upper left',fontsize=10)
    axes[3].plot(time,c1,'.', label='Total overlap (p1 overlaps p2)')
    axes[3].plot(time,c2,'.',label='Total overlap (p2 overlaps p1)')
    axes[3].plot(time,c3,'.',label='Partial overlap (p1 inner, p2 outer)')
    axes[3].plot(time,c4,'.',label='Partial overlap (p2 inner, p1 outer)')
    axes[3].legend(loc='upper left',fontsize=10)
    axes[3].set_ylim([0,1.2*max(max(c1),max(c2),max(c3),max(c4))])

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
