# Defect-Device-Degradation (DDD or "D3") Model, electrical degradation module
# Simulates solar cell power output degradation based on externally-calculated sodium migration profiles in silicon
# 01/15/2020
# ieng6/na299x/na299x/DB/Guillaume/Solar_cell_Al_BSF/Al_BSF

#__author__ = "Guillaume von Gastrow"
#__copyright__ = "Copyright 2020, SolEIL"
#__license__ = "MIT"
#__version__ = "1.0"
#__email__ = "guillaume.von.gastrow@alumni.aalto.fi"

# PIDmodel.py
# Module containing PID model functions

import os
import re # Regular Expression package
import pdb # debugger
import numpy as np
import matplotlib.pyplot as plt # plotting package
#import readh5conc
import csv
import h5py
from scipy import signal
from scipy import interpolate
import matplotlib.animation as manimation
import math

import runSent
import DatAnalysis

### batchshunt: Function running each sde and sdevice file consecutively and returning the IV curves.
### This functions runs the electrical module of the DDD model
### Takes as argument the conductivity profile as a function of time (time,sigma), and a Sentaurus simulation is performed for each time point
# Depth of the conductivity profile as a function of time (assumes shunt of fixed conductivity increasing in depth in the stacking fault)
# for each step, the function will create a sde file with updated shunt depth, the corresponding sdevice file, and run them
# Arguments:
# - batchdir:		Name of the directory where the generated sde and sdevice files will be saved. Example: "./test_dir"
# - changedir:		If the directory already exists, 1 will create a new one while 0 will save in the same one. By default changedir=1.
# - Temp:			Temperature (C)
# - mseg:			Value of the segregation coefficient of sodium from Si bulk to the stacking fault
# - clathrate_file	Name of the file containing the fit of clathrate resistivity as a function of Na to Si ratio
# - h5file:			Name of the h5py file containing the sodium migration profiles at each time point. eg "FOR_newNaprofile.h5"
# - folderpath:		Path to the folder containing Sentaurus templates files and where all data will be saved
# - startstep:		Number of the step where simulations will start. Should be 0 unless a different starting step is wanted (for instance if previous steps have already been run before).
# - endstep:		Number of the step where the simulations will end. If set to 0, simulations will run at all time points in the h5file containing the sodium migration profiles
# - skipNB:			Number of time points to skip in the sodium migration dataset at each iteration (0 by default). If set to 1, every second point will be run. If set to 2, every third point. etc.
# - sdetemplate		Name of the sdevice template file without .cmd extension (e.g. "sde_dvs")
# - sdevicetemplate	Name of the sdevice template file without .cmd extension (e.g. "sdevice_light_des")

#Naprofilename="single_layer_D1=4E-16cm2ps_D2=1E-15cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-04_m1.0e+00_pnp.h5"
Naprofilename="two_layers_D1=4E-16cm2ps_D2=1E-14cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-10_m1.0e+00_pnp.h5"
def runPIDsim(batchdir=".//testdir",changedir=1, Temp=60, mseg=10, use_SRV=False, clathrate_file="clathrates_cond.csv", h5file=Naprofilename, folderpath="/home/linux/ieng6/na299x/na299x/DB/Guillaume/Solar_cell_Al_BSF", startstep=0, endstep=0, skipNB=0, sdetemplate="sde_dvs", sdevicetemplate="sdevice_light_des"):
	# if defining end step, note that the simulation will stop one step before (for instance will stop at 5 if endstep=6)
	# "FOR_newNaprofile.h5"
	# "FOR_JV_85C.h5"	
	skip=skipNB+1 # skip is used in the division of the time point number (if skip=1, all points will be run)
	
	# Manage directories for file saving
	if changedir: # if true, a new directory will be created if the directory name already exists. If false, files will be saved in the existing directory
		if not os.path.exists(batchdir): 	# if the directory does not exist, create it
			os.makedirs(batchdir)
		else: 								# if the directory exists, change the name
			p=2
			while os.path.exists(batchdir): # correct this as if won't work beyond 9
				if batchdir[-1].isnumeric():
					if int(batchdir[-1])==p: # if the directory name already exists with the same number, increase the number (eg testdir3 becomes testdir4) (won't work beyond 9, can solve that by using regexp function to cut at the underscore and compare the numbers)
						p=p+1
						batchdir=batchdir[0:-1]+str(p)
						print("creating new directory name with different number")
				else:
					batchdir=batchdir+str(p) # if the director name exists but without a number, just add the number at the end
					print("creating new directory name")
			os.makedirs(batchdir)
	else:
		strdir='Using the existing directory '+batchdir
		print(strdir)
	
		#pdb.set_trace()

	# Name of the shunt to replace
	dshuntname="dshunt1"
		
	# #Open h5 file (file FOR_JV_85C.h5)
	# hf		= h5py.File(h5file, 'r')
	# ct		= hf.get('si_ct')
	# x		= hf.get('si_x')
	# time 	= hf.get('time')
	
	# for t in enumerate(time):
		# ct_field='ct_'+str(t[0]*2)
		
	# c at time i=0
	# c0 = ct[:,0]
	
	#Open h5 file (file for full stack file)
	hf		= h5py.File(h5file, 'r')
	time 	= hf.get('time')
	ct		= hf.get('L2/concentration') # L2 is the Si layer, L1 is the SiNx layer
	# ct_10=hf.get('si/concentration/ct_10')
	x		= hf.get('L2/x')
	
	#Open h5 file (file FOR_newNaprofile.h5)
	##hf		= h5py.File(h5file, 'r')
	##time 	= hf.get('time')
	##ct		= hf.get('si/concentration')
	# ct_10=hf.get('si/concentration/ct_10')
	##x		= hf.get('si/x')
	
	# c at time i=0
	c0=ct['ct_0'][:]
	
	# find nb of time steps from the h5 file
	# nbsteps=ct.shape[1]
	nbsteps=time.shape[0]
	
	# set last time step of the simulation in case it is not set by the user in the function arguments
	if endstep==0 or endstep<startstep:
		endstep=nbsteps
		print("Last simulation step set to the total nb of steps")
		
	# Create a log file in the simulation directory
	DatAnalysis.createlog(batchdir,Temp,mseg,clathrate_file, h5file, startstep, endstep, skipNB, sdetemplate, sdevicetemplate)
	
	#factor=50 # factor = 100 to have a final depth of 600 nm, as in "For_JV_85C.h5" the final depth is 6 nm. NOTE: Later, with the correct Na profiles, it should be just one.
	factor=1
	xdiff=(x-x[0])*factor # Difference because x does not start at 0
	#pdb.set_trace()
	
	# for each time point
	for i in range(startstep,endstep):
				# cond=condmodel(ct[:,i], Temp, clathrate_file, mseg)
		ct_field='ct_'+str(i) # field corresponding to profile at time i
		
		if np.mod(i,skip)==0 or endstep-startstep==1: # skip is related to the number of points to skip, and is 1 by default. If there is only 1 step (endstep-startstep==1), run it in any case.
		# if np.mod(i,2)==0 or endstep-startstep==1: # sample datapoints by only running even iterations (t0, t2, t4, t6, etc). If there is only 1 step (endstep-startstep==1), run it in any case.
			
			# Calculate conductivity profile in the shunt based on the Na concentration profile.
			cond=condmodel(ct[ct_field][:], Temp, mseg, clathrate_file)
			
			# Save conductivity profile in a .plx file
			condfilename="/conductivity_t_"+str(int(time[i]))+".plx"
			condfilepath=batchdir+condfilename # name of the conductivity file at this time point
			fp=open(condfilepath,'w')
			fp.write("\"MetalConductivity\"\n")
			for xloop,condloop in zip(xdiff,cond):
				line=str(xloop)+"\t"+str(condloop)+"\n" 
				fp.write(line)
				#pdb.set_trace()
			print("File "+condfilename+" created.\n")
			fp.close()
			#pdb.set_trace() #check that the conductivity file has been created
			# create new sentaurus files with modified shunt depth and conductivity
			
			## Create a new parameter file with updated SRV and save it in the simulation directory
			newParamfile=SRVparam(batchdir, ct[ct_field][:], Temp, time[i], use_SRV)
			
			# modify the sde and sdevice files and save them
			# NOTE: Includes shunt modification in the sde command file and modification of the parameter file name in the sdevice command file
			[newSDEname,newSDEVICEname]=changefiles(sdetemplate,sdevicetemplate,condfilename,newParamfile,ct[ct_field][:],xdiff,time[i],batchdir,dshuntname,mseg)
			
			# run Sentaurus
			#pdb.set_trace()
			
			#run Sentaurus Device Editor
			runSent.run_sde(newSDEname)
			
			# Check if sde execution caused an error, and if so stop
			if(errorcheck()):
				raise NameError('\n*******************\nError during sde execution. Stopping program.\n******************')
			
			# run sdevice
			runSent.run_sdevice(newSDEVICEname)
		
		# End of the loop for each time point
		
	# Close h5 file
	hf.close()
	
	finalstring="\nData in "+batchdir+".\n"
	print(finalstring)
	
	
	## Implementation of the conductivity model based on the interpolation of clathrate data.
	# Parameters:
	# - c is the sodium concentration in the Si bulk
	# - mseg is the segregation coefficient of sodium into stacking faults
	## Model simplifications:
	#- The Na to Si ratio in the stacking fault is obtained from the ratio between Na concentration and Si concentration in the bulk of a perfect crystal (does not consider the specific geometry of a stacking fault)
	#- Conductivity is calculated based on depth-resolved Hall-effect measurements of mobility and carrier density in Na-implanted Si (Korol et al)
def condmodel(c, Temp, mseg, clath_file="j"):

	cSi=5e22 # atomic density in Si (cm-3)
	#pdb.set_trace()
	# Import file containing fitting parameters for clathrate conductivity at different temperatures
	#with open(clath_file,newline='') as csvfile:
		##readline=csv.reader(csvfile,delimiter=',', quotechar='|',quoting=csv.QUOTE_NONNUMERIC)
	#	readline=csv.reader(csvfile)
		##for row in readline:
	#	pdb.set_trace()
		
	cshunt=c*mseg # Na concentration in the shunt
	
	## Clathrate model, not used because not realistic at our Na concentrations
	if False: # Skip this section with a trivial false condition
		Na_Si_ratio=cshunt/cSi # Model for Na density in shunt: ratio between Na density in shunt and Si bulk density in a perfect crystal
		# Add condition [Na]/[Si]=1 if higher than Si density
		
		# Later import the coefficients directly from the conductivity .csv file
		# coefficients at 60 deg C
		#a=24.960488582173806
		#b=-1.7580985174592794
		
		# coefs 70 C
		a=24.496024405605336
		b=-1.681009278808291
		
		# Calculate conductivity profile
		sigma=10**(a*Na_Si_ratio+b) # S/cm
		
		
	## Model based on implantation data
	# Korol, V. M. "Sodium ion implantation into silicon." Physica status solidi (a) 110.1 (1988): 9-34.
	# Fitting of coefficients in Extract_NaImp.py
	coord=-11.144769029961262
	slope=0.717839509854622
	
	sigma=(10**coord)*cshunt**slope # S/cm
	
	return sigma
	
	

# Function modifying the Sentaurus files to include updated shunt depth and modified external files (including the .par file)
# Creates a shunt with spatially varying conductivity based on doping profiles calculated at each time
# Inputs: 
# - sdetemplate:		Name of the Sentaurus Structure Editor file template without the .cmd extension (must be in the current folder)
# - sdevicetemplate:	Name of the Sentaurus Device file template without the .cmd extension (must be in the current folder)
# - condfilename:		Name of the updated .plx conductivity file
# - newParamfile:		Name of the updated .par parameter file
# - c:					Sodium concentration from the h5 file (list)
# - x:					depth from the h5 file (list)
# - time:				time from the h5 file (int or float)
# - dshuntname: 		name of the shunt to modify, as defined in the sde template
# - newfolderpath:		path of the directory where the newly generated Sentaurus files will be saved
# Example: runSentaurus.changefiles("sde_dvs","sdevice_light_des",2,"./test_dir")
def changefiles(sdetemplate, sdevicetemplate, condfilename, newParamfile, c, x, time, newfolderpath, dshuntname, mseg):
	# new name for the mesh file
	meshname="n_t"+str(int(time))
	
	# Calculate segregation coefficient at each depth and the depth of the shunt
	cseg=c*mseg
	L=len(cseg)
	#pdb.set_trace()
	k=0
	
	### ATTENTION: need to modify this as the Na concentration is not a meaningful value here. The shunt depth is highly sensitive to the conductivity model and to mseg.
	while(cseg[k]>0.1 and k<L-1): # Consider only the part of the [Na] profile in the stacking fault that is at least 1 cm-3 (to remove useless extra points). TO CHANGE.
		k=k+1
	
	##########################################################################""
	#### ATTENTION! This should be much more reliable than using cseg[k]>0.1, to be tested.
	# Or maybe we shouldn't limit the shunt depth at all.
	#while(PIDmodel.conmodel(cseg[k],60,mseg) > 1e-3 and k<L-1): # Limit the shunt depth to conductivities higher than 1e-3 S/cm (assuming 1e-3 is lower than the limit conductivity for shunting)
		#k=k+1
	###########################################################################"
		
	shdepth=x[k]-x[0] # Depth of the profile in um. The depth does not start at 0 so subtract x[0]
	
	# Define an arbitrary shunt depth in case of very low concentrations
	# that would lead to a shunt depth of 0 um and cause Sentaurus to crash
	if(shdepth==0):
		shdepth=0.1
		
	## NOTE: CHECK IF IT MATTERS IF THE EXTERNAL PROFILE IS DEEPER THAN THE DEPTH DEFINED HERE. CURRENTLY IT IS.
	#pdb.set_trace()
	
	########### open .cmd files to find number of lines ###########
		
	# open sde file to find number of lines
	fp=open(sdetemplate+".cmd",'r')
	count_sde = len(fp.readlines(  ))
	fp.close
	
	fp=open(sdevicetemplate+".cmd",'r')
	count_sdevice = len(fp.readlines(  ))
	fp.close
	
	linelist_sde=[0]*count_sde # preallocate list memory for original cmd file
	linelist_sdevice=[0]*count_sdevice # preallocate list memory for original cmd file
	newlinelist_sde=[0]*count_sde # preallocate list memory for cmd file with modified shunt depth
	newlinelist_sdevice=[0]*count_sdevice
	
	########### modify the sde file ############
	
	fp=open(sdetemplate+".cmd",'r') # reopen to return to beginning of the file
	k=0
	sdesuccesscount=0
	# save each line into a list
	while 1:
		dataline=fp.readline()
		if dataline=="": # check if the string is empty, meaning end of file
			break
		linelist_sde[k]=dataline # create a list containing one line at each index
		newlinelist_sde[k]=dataline # copy the file into a new list array
		
		# check where the shunt depth is defined and replace with updated shunt depth
		if re.search("define " + dshuntname, linelist_sde[k]):
			print("shunt line found")
			newlinelist_sde[k]="(define "+ dshuntname + " " + str(shdepth) +")\n" # replace by new depth value
			
		# Check where the conductivity file is defined in sde file and replace the path with the correct one
		#(define ShuntCondFile "./conductivity_test.plx")
		if re.search("define ShuntCondFile",linelist_sde[k]):
			print("External shunt file line found")
			newlinelist_sde[k]="(define ShuntCondFile " + "\"" + newfolderpath + condfilename +"\")\n" # replace by new depth value
			
		# find line where mesh and ouptut file names are defined and change their names at each simulation time
		if re.search("\"nodnum1\"",linelist_sde[k]): # if the mesh definition is found
			print("mesh line found")
			newlinelist_sde[k]=re.sub("nodnum1",newfolderpath+"//"+meshname,linelist_sde[k]) # the mesh file is saved in the new folder created for the current batch simulation
			sdesuccesscount=sdesuccesscount+1
			#pdb.set_trace()
			
		k=k+1
	fp.close() #close sdevice template file
	print("SDE file loaded")
	
	########### modify the sdevice file #############
	
	fp=open(sdevicetemplate+".cmd",'r')
	k=0
	successdevice=0
	# save each line into a list
	while 1:
		dataline=fp.readline()
		if dataline=="": # check if the string is empty, meaning end of file
			break
		linelist_sdevice[k]=dataline # create of list containing one line at each index
		newlinelist_sdevice[k]=dataline # copy the file into a new list array
		
		## Update the node names (called "nodnum1" in template)
		if re.search("nodnum1",linelist_sdevice[k]): # if the mesh definition is found
			print("mesh line found")
			newnodename="n_t"+str(int(time))
			newlinelist_sdevice[k]=re.sub("nodnum1",newfolderpath+"//"+newnodename,linelist_sdevice[k]) # replace the nodenumber "nodnum1" in the template by n_t1, n_t2, etc.
			successdevice=successdevice+1
			#pdb.set_trace()
			
		## Update the name of the parameter file (called "sdevice" in the template)
		if re.search("sdevice",linelist_sdevice[k]): # if the parameter file is found
			print("parameter file line found")
			newlinelist_sdevice[k]=re.sub("sdevice.par",newParamfile,linelist_sdevice[k]) # replace the nodenumber "nodnum1" in the template by n_t1, n_t2, etc.
			successdevice=successdevice+1
		
		k=k+1
	fp.close()
	
	################# save the files #################
	
	# save the files under a different name (_t0, _t1, _t2, etc)
	# file path
	sdefilepath=os.path.join(newfolderpath,sdetemplate)
	sdevicefilepath=os.path.join(newfolderpath,sdevicetemplate)
	
	newSDEname=sdefilepath+"_t"+str(int(time))+".cmd"
	f=open(newSDEname,"w+")
	for i in range(len(newlinelist_sde)):
		f.write(newlinelist_sde[i])
	f.close()
	print("File "+newSDEname+" created.")
	
	newSDEVICEname=sdevicefilepath+"_t"+str(int(time))+".cmd"
	f=open(newSDEVICEname,"w+")
	for i in range(len(newlinelist_sdevice)):
		f.write(newlinelist_sdevice[i])
	f.close()
	
	print("File "+newSDEVICEname+" created.")
	#pdb.set_trace()
	
	# Check whether all the names where updated in sde and sdevice files
	if sdesuccesscount<1: # in the template, only 1 line should contain "nodnum1", the line mesh is built (only the number of lines are counted)
		warning_message="\nATTENTION:\nMesh naming may have failed (check generated sde file at time" + "t"+str(int(time)) +").\nWould you like to continue?\ny/n\n"
		choice=''
		while choice !='y' and choice !='n':
			choice=input(warning_message)
		if choice=='n':
			raise NameError('\n\n*******************************\nSimulations cancelled by user.\n**************************\n')			
	
	if successdevice<5: # the sde file uses 5 files named at each simulation time: the mesh file from sde, the optical generation, the output data file (.plt, output), the results file (.tdr, output) and the parameter file (sdevice.par)
		warning_message="\ATTENTION:\nFile naming may have failed (check generated sdevice file at time" + "t"+str(int(time)) +").\nWould you like to continue?\ny/n\n"
		choice=''
		while choice !='y' and choice !='n':
			choice=input(warning_message)
		if choice=='n':
			raise NameError('\n\n*******************************\nSimulations cancelled by user.\n**************************\n')			
	
	return newSDEname, newSDEVICEname
	
# Function giving the parameterization of surface recombination velocity as a function of the surface Na concentration
def SRVparam(batchdir, cNa, T, time, use_SRV):
# use_SRV	Boolean. if use_SRV=False, S0 is set to zero
# batchdir	Directory where the simulation files are saved
# cNa		List containing Na concentration as a function of depth (cm-3)
# T			Temperature (C) (int or float)
# time		Time at which the sodium profile was calculated (int or float)

	if not use_SRV:
		S0=0 # cm/s
	else:
		# Calculate surface recombination velocity
		# Fitting values according to the phosphorus parameterization by Altermatt et al, Journal of App Phys 92, 3187, 2002.
		S1=500 # cm/s
		S2=60 #cm/s
		gamma1=0.6
		gamma2=3
		N1=1e10 # cm-3 (modified from Altermatt et al)
		N2=1e10 # cm-3 (modified from Altermatt et al)
		
		# Parameterization of the surface recombination velocity
		S0=S1*(cNa[0]/N1)**gamma1+S2*(cNa[0]/N2)**gamma2 # Altermatt et al, Journal of App Phys 92, 3187, 2002
		
		# Limit S0 to the thermal velocity of electrons
		
		me=9.1e-31 # electron mass, kg
		kB=1.38e-23 # J.K-1
		TK=T+273.15
		
		vth=100*math.sqrt(3*kB*TK/me) # (cm/s) thermal energy for non-relastivistic electrons E=df*kB*T where df number of degrees of freedom, and E=m*v^2
		
		if S0>vth: # cm/s
			S0=vth # cm/s
	
	# Limit to 5 significant digits
	S0_5g=format(S0, '1.4e')
	
	# Modify the sdevice.par template file to include this surface recombination velocity
	# and save the modified .par file in batchdir
	newline=S0_5g+" ,\t"+S0_5g
	newParamfile=DatAnalysis.replace_line("sdevice","par","S0_val , S0_val",newline,batchdir,time,"SRV value")
	
	return newParamfile

# check errfile to see whether sde execution caused an error
# Returns
# True if there is an error
# False if there is no error
def errorcheck():
	with open('errfile.txt', 'r') as f:
		err_flag=f.read()
	#pdb.set_trace()
	return bool(int(err_flag))
		
	
#######################################################
### Functions below were for testing purposes only

	# Function giving depth of the shunt as a function of time
def shuntcond(nbsteps):
	sigma=np.zeros((nbsteps,2))
	sigma[0,1]=0.2 # Depth of the shunt in um at t=0 (correct this part if need to redo data analysis, before was 0.8)
	for i in range(1,nbsteps): #start loop at t=1
		sigma[i,0]=i # time in arbitraty units as 0,1,2,3,...,nbsteps-1
		sigma[i,1]=sigma[i-1,1]+0.1
		# here add a condition to reduce the conductivity in the parameter file if the pn junction has been reached
	return sigma
	
	# Function to plot conductivity curves from all .plx files found in a folder
def plotcond(folderpath):
	print("ok")