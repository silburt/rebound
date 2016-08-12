import matplotlib.pyplot as plt
import numpy as np

fos = open('energy.txt', 'r')
data = np.loadtxt(fos, delimiter=',')
fig, a = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(10,8))

a[0].plot(data[:,0],data[:,1],'.')
a[0].set_ylabel('energy')
a[0].set_yscale('log')
a[1].plot(data[:,0],data[:,2],'.',label='planet 1')
a[1].plot(data[:,0],data[:,3],'.',label='planet 2')
a[1].set_ylabel('a')
a[2].plot(data[:,0],data[:,4],'.')
a[2].set_ylabel('N')
a[3].plot(data[:,0],data[:,6],'.')
a[3].set_ylabel('HSF')

a[0].set_xscale('log')

plt.show()