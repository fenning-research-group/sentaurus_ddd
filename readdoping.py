# reads h5py doping profile

import matplotlib.pyplot as plt
import pdb
import h5py

filename="FOR_newNaprofile.h5"

hf=h5py.File(filename,'r')

time 	= hf.get('time')
ct		= hf.get('si/concentration') # si
# ct_10=hf.get('si/concentration/ct_10')
x		= hf.get('si/x')


#pdb.set_trace()
# c at time i=0
c21=ct['ct_21'][:] #time at 34560 s where it seems that PID starts

plt.ion()
plt.figure()

xreal=x-x[0]

plt.plot(xreal,c21)
plt.yscale('log')
plt.ylim(1e5,1e20)
plt.xlim(0,0.8)
plt.xlabel('Depth (um)')
plt.ylabel('Sodium concentration (cm^-3)')
plt.subplots_adjust(left=0.15,top=0.9,wspace=0.4) # adjust space around plot

hf.close()