# Plot simulation results with and without SRH

import matplotlib.pyplot as plt
import pdb
import PIDmodel
import numpy as np
import h5py

"""
# Compare simulation results with and without SRH
[res_noSRH,timepts_noSRH,effpts_noSRH,timenew_noSRH,effsmooth_noSRH]=PIDmodel.DatAnalysis.batchanalysis('./20200127_curve_NO-SRH2_saved',startstep=30,endstep=80)
[res_SRH,timepts_SRH,effpts_SRH,timenew_SRH,effsmooth_SRH]=PIDmodel.DatAnalysis.batchanalysis('./20200127_curve_SRH2_saved',startstep=30,endstep=80)

plt.ion()
plt.figure()

# xp points
plt.plot(timepts_noSRH/3600,effpts_noSRH,marker='s',linestyle='',color='k',fillstyle='none',markersize=12,label='noSRV')
plt.plot(timepts_SRH/3600,effpts_SRH,marker='^',linestyle='',color='k',fillstyle='none',markersize=12,label='SRV')

# interpolated curve
plt.plot(timenew_noSRH/3600,effsmooth_noSRH,'--b',linewidth=1.3,dashes=(5,6)) #plot interpolated curve
plt.plot(timenew_SRH/3600,effsmooth_SRH,'--b',linewidth=1.3,dashes=(5,6)) #plot interpolated curve

plt.subplots_adjust(left=0.15,top=1,wspace=0.4) # adjust space around plot

plt.ylabel('Efficiency (%)')
plt.xlabel('Migration time (h)')

ax=plt.gca()
leg=ax.legend()

pdb.set_trace()
"""
"""
## plot phosphorus/boron doping profile
import scipy.special as sci
import numpy as np
x=np.linspace(0,1,500)
x=np.linspace(0,1.5,500)
c=1e19*(1-sci.erf(x/0.2))
import matplotlib.pyplot as plt
plt.ion()
plt.figure()
plt.plot(x,c)
plt.yscale('log')
plt.ylim(1e15,1e20)
plt.xlim(0,0.8)
plt.ylabel('Emitter doping concentration (cm$^{-3}$)')
plt.xlabel('Depth ($\mu$m')
plt.xlabel('Depth ($\mu$m)')
plt.subplots_adjust(left=0.15,top=0.9,wspace=0.4) # adjust space around plot

ax=plt.gca()

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
	ax.get_xticklabels() + ax.get_yticklabels()):
	item.set_fontsize(19)

plt.subplots_adjust(left=0.15,top=0.95,bottom=0.15,wspace=0.4) # adjust space around plot

"""

## Plot sodium profiles at different times
#Open h5 file (file for full stack file)
#Naprofilename="single_layer_D1=4E-16cm2ps_D2=1E-15cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-04_m1.0e+00_pnp.h5"
Naprofilename="two_layers_D1=4E-16cm2ps_D2=1E-14cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-10_m1.0e+00_pnp.h5"
hf		= h5py.File(Naprofilename, 'r')
time 	= hf.get('time')
ct		= hf.get('L2/concentration') # L2 is the Si layer, L1 is the SiNx layer
	# ct_10=hf.get('si/concentration/ct_10')
x		= hf.get('L2/x')
vfb=hf.get('vfb')
plt.ion()
plt.figure()
for ii in range(len(ct)):
	# cond=condmodel(ct[:,i], Temp, clathrate_file, mseg)
	ct_field='ct_'+str(ii) # field corresponding to profile at time i
	plt.semilogy(x,ct[ct_field])
	
plt.ylim(1e0,1e18)
plt.xlim(0.07,0.8)

plt.xlabel('Depth (um')
plt.ylabel('Sodium concentration (cm$^{-3}$)')

plt.subplots_adjust(left=0.2,top=0.95,bottom=0.15,wspace=0.4) # adjust space around plot

ax=plt.gca()

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
	ax.get_xticklabels() + ax.get_yticklabels()):
	item.set_fontsize(19)
	
#plt Vfb as a function of time
plt.figure()
plt.plot(time[()]/3600,vfb)

#plt.ylim(1e0,1e18)
#plt.xlim(0.07,0.8)

plt.xlabel('Time (hrs)')
plt.ylabel('Vfb shift (V)')

plt.subplots_adjust(left=0.2,top=0.95,bottom=0.15,wspace=0.4) # adjust space around plot

ax=plt.gca()

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
	ax.get_xticklabels() + ax.get_yticklabels()):
	item.set_fontsize(19)

pdb.set_trace()
hf.close()

if False:
	## Plot sodium profile at the threshold for PID in the SRH/noSRH comparion from 2020/01/27
	Naprofilename="single_layer_D1=4E-16cm2ps_D2=1E-15cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-04_m1.0e+00_pnp.h5"

	ct_field='ct_68' # point corresponding to 32640 s, the closest to 9 hours
	hf		= h5py.File(Naprofilename, 'r')
	time 	= hf.get('time')
	ct		= hf.get('L2/concentration') # L2 is the Si layer, L1 is the SiNx layer
		# ct_10=hf.get('si/concentration/ct_10')
	x		= hf.get('L2/x')

	plt.semilogy(x,ct[ct_field],color=(0.235,0.643,0.45),linewidth=2)
	# RGB 60,164,115
	plt.ylim(1e5,1e18)
	plt.xlim(0.07,0.8)

	plt.xlabel('Depth (um')
	plt.ylabel('Sodium concentration (cm$^{-3}$)')

	plt.subplots_adjust(left=0.2,top=0.95,bottom=0.15,wspace=0.4) # adjust space around plot

	ax=plt.gca()

	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
		ax.get_xticklabels() + ax.get_yticklabels()):
		item.set_fontsize(19)
