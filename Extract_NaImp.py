# code to extract and plot conductivity data from Korol, V. M. "Sodium ion implantation into silicon." Physica status solidi (a) 110.1 (1988): 9-34.
# Data from curve 1 in figure 13

import matplotlib.pyplot as plt
import pdb
import numpy as np
import csv
import os
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
import statistics as stat

import PIDmodel

mseg=100
dSi=5e22 # cm-3, approx. atomic density is silicon

def csvextract(csvpath):
# Function extracting data from a 2-column, comma separated .csv file
# and returning interpolated values

	#initialize list
	xdata=[] # depth (nm)
	ydata=[] # carrier density (cm-3)

	with open(csvpath,newline='') as cf:
		reader=csv.reader(cf,delimiter=',', quotechar='|')
		for row in reader:
			xdata.append(float(row[0]))
			ydata.append(float(row[1]))
			
	return xdata,ydata
	
def csvinterp(csvpath1,csvpath2):
	# Function interpolating csv data output from csvextract with same number of points for conductivity and mobility

	# Extract data from .csv file and interpolate
	[xdata1,ydata1]=csvextract(csvpath1)
	[xdata2,ydata2]=csvextract(csvpath2)
	
	# Find left boundary for interpolation (largest of the minima)
	xmin1=min(xdata1)
	xmin2=min(xdata2)
	xinterpmin=max(xmin1,xmin2)
	
	# Find right boundary for interpolation (smallest of the maxima)
	xmax1=max(xdata1)
	xmax2=max(xdata2)
	xinterpmax=min(xmax1,xmax2)

	# Array used for interpolation
	xinterp=np.linspace(xinterpmin,xinterpmax,200) # later use min and max of both curves

	## interpolate curves with the same x-array to have the desired number of points
	y_interpobj1=interp1d(xdata1,ydata1)
	yinterp1=y_interpobj1(xinterp)
	
	y_interpobj2=interp1d(xdata2,ydata2)
	yinterp2=y_interpobj2(xinterp)
	
	return xdata1,ydata1,xdata2,ydata2,xinterp,yinterp1,yinterp2
	
# Function to split a curve so that we only have one ordinate per point
def splitcurve(x,y):
	k=0
	while x[k+1]>x[k]:
		k=k+1
	x1=x[0:k]
	x2=x[k+1:]
	y1=y[0:k]
	y2=y[k+1:]
	
	return x1,y1,x2,y2
	
def fun_exp(x,cNa,sigma):
	return (sigma-x[0]*cNa**x[1])
	
def sigmafunc(CNa,coord,slope):
	return 10**coord*CNa**slope

## Path to the experimental data extracted from Korol et al, fig 13, Na implantation in Si.
path='//home//linux//ieng6//na299x//na299x//DB//Guillaume//cell_AlBSF//Al_BSF//Na_implantation_data'

mobilityfile_c1='Korov_curve1_mobility.csv' # mobility, curve 1 (cold implantation)
carrierfile_c1='Korov_fig13_curve1_n.csv' # carrier density, curve 1 (cold implantation)

mobilityfile_c2='Korov_fig13_curve2_mobility_hot.csv'
carrierfile_c2='Korov_fig13_curve2_n_hot.csv'

mupath_c1=os.path.join(path,mobilityfile_c1)
carrierpath_c1=os.path.join(path,carrierfile_c1)

mupath_c2=os.path.join(path,mobilityfile_c2)
carrierpath_c2=os.path.join(path,carrierfile_c2)

depth_mu=[] # depth (nm)
depth_n=[] # depth (nm)
carrierdens=[] # carrier density (cm-3)
mobility=[] # mobility (cm2V-1s-1)

#constants
q=1.60219e-19 # C

# Extract data from .csv file and interpolate
[depth_mu_rawc1,mu_rawc1,depth_n_rawc1,n_rawc1,interpdepth_c1,mu_interp_c1,n_interp_c1]=csvinterp(mupath_c1,carrierpath_c1)
[depth_mu_rawc2,mu_rawc2,depth_n_rawc2,n_rawc2,interpdepth_c2,mu_interp_c2,n_interp_c2]=csvinterp(mupath_c2,carrierpath_c2)

# Plot raw experimental data
plt.ion()
plt.figure()
plt.semilogy(depth_n_rawc1,n_rawc1,'+')
plt.semilogy(depth_mu_rawc1,mu_rawc1,'+')

plt.ylim(1,1e18)
plt.xlim(0,575)

# Plot interpolated data
plt.semilogy(interpdepth_c1,n_interp_c1)
plt.semilogy(interpdepth_c1,mu_interp_c1)
plt.xlabel('Depth (nm)')
plt.ylabel('Carrier density (cm-3) or mobility (cm2/V/s)')

# Conductivity (cm-1Ohm-1)
sigma_c1=n_interp_c1*mu_interp_c1*q # cm-1*Ohm-1
sigma_c2=n_interp_c2*mu_interp_c2*q

plt.figure()
plt.loglog(n_interp_c1/mseg,sigma_c1,'+',label='Implantation curve1 (lowT)') # Plot conductivity in the shunt as a function of [Na] in the Si bulk (n_interp_c1/mseg).
plt.loglog(n_interp_c2/mseg,sigma_c2,'+',label='Implantation curve2 (highT)') # Plot conductivity in the shunt as a function of [Na] in the Si bulk (n_interp_c1/mseg).
plt.xlabel('Na concentration (cm-3)')
plt.ylabel('conductivity (S/cm)')
plt.ylim(0.5,100)
plt.xlim(1e15,3e18)

## Plot clathrate conductivity model as a function of sodium concentration
bool=True
if bool: # if bool=false, ignore this section
	c=np.logspace(10,21,1000)
	plt.loglog(c,PIDmodel.condmodel(c,85,mseg),label='clathrate') # segregation coeff of 120
	#plt.loglog(c,PIDmodel.condmodel(c,85,mseg*100),label='mseg=1e4') # segregation coeff of 120
	#plt.loglog(c,PIDmodel.condmodel(c,85,mseg*1e3),label='mseg=1e5') # segregation coeff of 120
	plt.xlim(1e10,1e22)
	plt.ylim(1e-3,1e5)
	plt.xlabel('Sodium Concentration in Si bulk (cm$^{-3}$)')
	plt.ylabel('Conductivity of the shunt (mseg=100) (S/cm)')

	ax=plt.gca()
	leg=ax.legend()

	plt.subplots_adjust(left=0.2,top=0.95,bottom=0.18,wspace=0.4) # adjust space around plot

	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
		ax.get_xticklabels() + ax.get_yticklabels()):
		item.set_fontsize(19)

	# Plot actual clathrate experimental points
	c1=3/136*dSi/mseg # corresponding Na concentration in the Si bulk
	p1=PIDmodel.condmodel(c1,85,mseg)
	plt.loglog(c1,p1,'+',label='Na3Si136')

	c2=8/46*dSi/mseg # corresponding Na concentration in the Si bulk
	p2=PIDmodel.condmodel(c2,85,mseg)
	plt.loglog(c2,p2,'+',label='Na8Si46')
	

# Separate each curve to have only a single ordinate per point
[n_interp_c1_1,sigma_c1_1,n_interp_c1_2,sigma_c1_2]=splitcurve(n_interp_c1[5:],sigma_c1[5:]) # start at 5 where the x values start increasing in the list
[n_interp_c2_1,sigma_c2_1,n_interp_c2_2,sigma_c2_2]=splitcurve(n_interp_c2,sigma_c2)

plt.figure()
plt.loglog(n_interp_c1_1,sigma_c1_1,'+',label='c1_1')
plt.loglog(n_interp_c1_2,sigma_c1_2,'+',label='c1_2')
plt.loglog(n_interp_c2_1,sigma_c2_1,'+',label='c2_1')
plt.loglog(n_interp_c2_2,sigma_c2_2,'+',label='c2_2')

#pdb.set_trace()

# Fit each curve with a function y=a*x^b. Then average to obtain final model.
# Exponential function to fit conductivity curves
x0=[1e-17,1]
xc1_1=least_squares(fun_exp, x0, args=(n_interp_c1_1, sigma_c1_1))
xc1_2=least_squares(fun_exp, x0, args=(n_interp_c1_2, sigma_c1_2))
xc2_1=least_squares(fun_exp, x0, args=(n_interp_c2_1, sigma_c2_1))
xc2_2=least_squares(fun_exp, x0, args=(n_interp_c2_2[30:], sigma_c2_2[30:]))

#pdb.set_trace()
yfitc1_1=xc1_1.x[0]*n_interp_c1_1**xc1_1.x[1]
#plt.loglog(n_interp_c1_1,yfitc1_1)

yfitc1_2=xc1_2.x[0]*n_interp_c1_2**xc1_2.x[1]
#plt.loglog(n_interp_c1_2,yfitc1_2)


yfitc2_1=xc2_1.x[0]*n_interp_c2_1**xc2_1.x[1]
#plt.loglog(n_interp_c2_1,yfitc2_1)

yfitc2_2=xc2_2.x[0]*n_interp_c2_2**xc2_2.x[1]
#plt.loglog(n_interp_c2_2[30:],yfitc2_2[30:],label='c2_2 fit')

plt.xlabel('Na concentration (cm$^{-3}$')
plt.ylabel('Conductivity (S/cm)')

# Average fits (non-weighted average for now, anyway it is a coarse fit)
coord_mean=stat.mean([np.log10(xc1_1.x[0]),np.log10(xc1_2.x[0]),np.log10(xc2_1.x[0]),np.log10(xc2_2.x[0])])
slope_mean=stat.mean([xc1_1.x[1],xc1_2.x[1],xc2_1.x[1],xc2_2.x[1]])

plt.loglog(n_interp_c2_2,10**coord_mean*n_interp_c2_2**slope_mean,label='global fitted sigma')

plt.legend()

plt.figure()
CNa=np.logspace(-0,22)
plt.loglog(CNa,sigmafunc(CNa,coord_mean,slope_mean))
plt.xlabel('Na concentration (cm$^{-3}$')
plt.ylabel('Average conductivity (S/cm)')

# Return the fitting coefficients of conductivity (S/cm) as a function of Na concentration (cm-3)
print(coord_mean);print(slope_mean)