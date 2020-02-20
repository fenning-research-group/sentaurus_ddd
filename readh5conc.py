import h5py
import pdb
import matplotlib.pyplot as plt
import numpy as np
from runSentaurus import condmodel

hf		= h5py.File('FOR_JV_85C.h5', 'r')
ct		= hf.get('si_ct')
x		= hf.get('si_x')
time 	= hf.get('time')

# c at time i=0
c0 = ct[:,0]

# c10=ct[:,10] # concentration at the 10th time point (first one is at t=0)
c3=ct[:,3] 

mseg=10
cseg=c3*mseg

k=0
while(cseg[k]>1): # Consider only the part of the [Na] profile in the stacking fault that is at least 1 cm-3
	k=k+1
	
depth=x[k]-x[0] # Depth of the profile in um. The depth does not start at 0 so subtract x[0]

timepts=ct.shape[1] # number of time points
#for i in range(timepts):
#	plt.plot(x[:]*1e3,ct[:,i])
#	plt.yscale('log')
#	plt.xlabel('Depth (nm)')
#	plt.ylabel('Na (cm-3)')
#	plt.ylim(1e15,1e22)
#	plt.xlim(70, 83)
#	
#plt.show()

#pdb.set_trace()

x0=(x[:]-x[0])*1e3
x0norm=x0/1e3

#test saving conductivity file
# # cond=condmodel(ct[:,10], 60)
		
		# # # Save conductivity profile in a text file
# # condfilename="./test_dir/conductivity_t_"+str(int(time[10]))+".plx" # name of the conductivity file at this time point
# # fp=open(condfilename,'w')
# # fp.write("\"MetalConductivity\"\n")
# # for xloop,condloop in zip(x0norm,cond): # loop over the values in x and cond
	# # line=str(xloop*100)+"\t"+str(condloop)+"\n" # *100 to have a final depth of 600 nm, as in "For_JV_85C.h5" the final depth is 6 nm
	# # fp.write(line)
	
# # print("File "+condfilename+" created.")
# # fp.close

hf.close()

pdb.set_trace()
sigma=condmodel(c10,50)

plt.plot(x0*100,sigma)
plt.xlabel("depth (nm)")
plt.ylabel("sigma (S/cm)")
plt.yscale('log')
plt.xlim(0, 700)
plt.ylim(1e-2, 1.5)

cratio=c10*mseg/5e22
plt.figure()
plt.plot(x0*100,cratio) # verify that the clathrate model is applied correctly
plt.xlabel("depth (nm)")
plt.ylabel("[Na]shunt/[Si]")
plt.yscale('log')
plt.xlim(0, 700)
plt.ylim(1e-8, 1.5)

plt.figure()
plt.plot(x0*100,c10) # verify that the clathrate model is applied correctly
plt.xlabel("depth (nm)")
plt.ylabel("[Na]")
plt.yscale('log')
plt.xlim(0, 700)
plt.ylim(1e15, 1e21)


plt.show()

pdb.set_trace()
