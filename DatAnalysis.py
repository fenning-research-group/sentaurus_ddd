# DatAnalysis.py
# Module containing data analysis functions used for the DDD model

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
from datetime import datetime


# function to plot the IV curve
# Arguments:
#		filename: name of the .plt file containing the data. Example: "n_t6600_light_des.plt"
#		plotcond=1 to plot, 0 otherwise
def analyzedata(filename,plotcond):

	#if condition=="dark":
	#	fid=open("nodnum1_dark_des.plt","r")
	#elif condition=="light":
	#	fid=open("nodnum1_light_des.plt","r")
	#else:
	#	raise Exception('wrong argument entered. Enter "dark" or "light".')	
	
	fid=open(filename,"r")
	rawdata="" #Create empty string
	while 1:
		dataline=fid.readline()
		if dataline=="": # check if the string is empty, meaning end of file
			break
		linestr=dataline.strip() # save dataline into a string and remove leading and trailing spaces
		rawdata=rawdata+" "+linestr # concatenate each line with the line string
	fid.close()
	
	# find the indices of the { brackets
	# xstart=re.search("{",data)
	# xstart.span()
	indstart=[m.span() for m in re.finditer("{",rawdata)] # finditer finds all the iterations of '{'
	
	# find the indices of the } brackets
	indend=[m.span() for m in re.finditer("}",rawdata)]
	
	# Make a dictionary 'data' containing the info section and the numerical data section as an array
	rawinfo=rawdata[indstart[0][0]+1:indend[0][0]-1] # names of datasets and of functions (includes extra spaces)
	rawvalues=rawdata[indstart[1][1]+1:indend[1][1]-1] # numerical values calculated by Sentaurus (includes extra spaces)
	
	# create the dictionary entries and remove spaces
	data=dict() # define dictionary
	data["info"]=rawinfo.strip() 
	data["values"]=rawvalues.strip()
	
	#pdb.set_trace()
	# Find indices of the [ and ] brackets in the info section
	# [ brackets
	sqindst=[m.span() for m in re.finditer("\[",rawinfo.strip())] # starting index of square bracket
	# ] brackets
	sqinden=[m.span() for m in re.finditer("\]",rawinfo.strip())] # ending index of square bracket
	
	# The dataset names are within the first brackets
	rawdatasets=data["info"][sqindst[0][0]+1:sqinden[0][0]-1]
	data["datasets"]=rawdatasets.strip()
	
	# Split the dataset field at the double quotes to find the number of output parameters (voltage, current, etc.)
	list_datasets=re.split("\" ",data["datasets"]) # split each time a double quote followed by a space in encountered
	nboutputs=len(list_datasets) #number of output parameters from sdevice
	
	# Split the value string list
	#list_values=re.split(" ",data["values"]) # separate with two spaces as this is the minimum spacing between consecutive values
	list_values=data["values"].split() # split without argument splits at white spaces
	
	# Convert the value string list into a numpy string array (vector)
	valvect_str=np.array(list_values)
	
	# convert numpy string array into a numpy float array
	valvect=valvect_str.astype(np.float)
	len_valvect=len(valvect)
	
	# Reshape the data array
	nblines=int(len_valvect/nboutputs)# nb of lines for the matrix
	#pdb.set_trace()
	valarray=np.reshape(valvect,(nblines,nboutputs))
	#pdb.set_trace()
	
	# Get current and voltage (column number depends on the defined outputs in the sdevice file)
	V=valarray[:,1] # voltage in col. index 1, ie 2nd column
	I=valarray[:,7] # voltage in col. index 7, ie 8th column
	
	if plotcond:
		plt.ion()
		plt.plot(V,I,'-+')
		plt.ylabel("Total current (mA/cm2)")
		plt.xlabel("Voltage (V)")
		plt.ylim(-50, 20)
		plt.show()
	
	return [I, V]
	
Naprofilename="single_layer_D1=4E-16cm2ps_D2=1E-15cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-04_m1.0e+00_pnp.h5"
#Naprofilename="two_layers_D1=4E-16cm2ps_D2=1E-14cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-10_m1.0e+00_pnp.h5"
def batchanalysis(directory='./2019-11-1_backup',sdevicetemplate="sdevice_light_des",h5File=Naprofilename,startstep=0,endstep=0):
	# analyzes data from several plt files
	# assumes that the .plt files are saved in directory
	#sdevicetemplate is the name of the Sentaurus device template file without the extension, for instance "sdevice_light_des"
	# example: runSentaurus.batchanalysis(6,"sdevice_light_des")
	
	# obtain depth of the shunt as a function of time
	#sigma=shuntcond(nbsteps)
	
	# Find time point corresponding to each time
	# Open h5 file
	#h5File="FOR_JV_85C.h5"
	hf= h5py.File(h5File, 'r')
	tme= hf.get('time')
	ti=tme[:]
	
	Lti=len(ti) # length of the time list
	
	# In case endstep is not defined by the user or is wrongly defined
	if endstep==0 or endstep<startstep:
		endstep=Lti
		print("Last simulation step set to the total nb of steps")
		
	#nbsteps=len(ti) # number of time steps
	nbsteps=endstep-startstep

	# Save h5py data
	time=[0]*nbsteps
	time[:]=ti[:]
	
	hf.close()
		
	# Define dictionary containing IV data. Each field contains a list.
	results=dict()
	# Create lists of NaN (so that 0's won't be plotted if some points are skipped)
	results["I"]=NaNlist(Lti)
	results["sdevicename"]=NaNlist(Lti)
	results["time"]=NaNlist(Lti)
	results["V"]=NaNlist(Lti)
	results["shdepth"]=NaNlist(Lti)
	results["Efficiency"]=NaNlist(Lti)
	
	#pdb.set_trace()
	# results["I"]=[0]*nbsteps # intialize the list to save current data
	# results["V"]=[0]*nbsteps
	# results["sdevicename"]=[0]*nbsteps
	# results["time"]=[0]*nbsteps
	# results["shdepth"]=[0]*nbsteps
	# results["Efficiency"]=[0]*nbsteps
	
	shdepth=[]# List containing depth of the shunt at each iteration
	
	plt.ion()
	plt.figure()
	for i in range(startstep,endstep):
		# Analyze data and store generated IV curve in a dictionary
		# The name of the file to analyze is obtained by adding the node name to the sdevice file name without the string "sdevice".
		# eg: sdevice_light_des.cmd becomes n_t1_light_des.plt
		
		#shdepth=sigma[i,1] # append shunt depth
		#timestamp="_t"+str(int(sigma[i,0]))
		
		basename=re.split("sdevice_",sdevicetemplate) # yields "light_des" from "sdevice_light_des", case of a template file "sdevice_light_des.cmd".
		plotfile=directory+"//"+"n_t"+str(int(time[i]))+"_"+basename[1]+".plt" # by default the .plt files containing the extracted data are saved with a name in the form "n_tk_basename.plt" where k is the time and basename is the string after "sdevice" in the sdevice file. Example: "n_t4_light_des.plt" based on the file "sdevice_light_des_tk.cmd".
		
		# Try openining the filename. Pass if it does not exist (in case some points were skipped in the simulations)
		plt_flag=1
		try:
			fid=open(plotfile,"r")
			fid.close()
		except Exception as e:
			print(e)
			warning_str='Could not open file '+plotfile+', it may not exist.\nSkipping the file in the analysis.\n'
			print(warning_str)
			plt_flag=0 # set the flag to 0 so the file won't be used in the analysis

		if plt_flag: # if the .plt file was found, extract data
			# Extract resulting IV curve
			print(plotfile)
			#pdb.set_trace()
			
			try: # try extracting data from the .plt file
				[I,V]=analyzedata(plotfile,0)
			except: # if the file is corrupted, skip this .plt file.
				# (happens for instance if the simulations do not converge)skip this .plt file.
				# The corresponding dictionary value will remain a NaN.
				print('*****File %s is corrupted, skipped in analysis\n' % plotfile)
			else: # if no exception is found, save the extracted data in the dictionary
				# Extract efficiency as a function of time step
				# save results in dictionary
				print('File %s is ok\n' % plotfile)
				results["I"][i]=I
				results["V"][i]=V
				results["sdevicename"][i]=plotfile
				#results["time"][i]=i
				results["time"][i]=time[i]
				results["Efficiency"][i]=findeff(V,I)
				#results["shdepth"][i]=shdepth
				
				# Plot the IV curves
				#plt.figure()
				plt.subplot(1,2,1)
				plt.plot(V,I)
			
	plt.ylabel("Total current (mA/cm$^2$)")
	plt.xlabel("Voltage (V)")
	plt.ylim(-35, 0)
	plt.xlim(0, 0.75)
	plt.rcParams.update({'font.size':16})
	
	# pdb.set_trace()
	
	# Convert lists to numpy arrays to apply a mask
	timearr=np.asarray(results["time"])
	effarr=np.asarray(results["Efficiency"])
	# Find where to mask NaN in the array
	effmask=np.isfinite(effarr)
	#pdb.set_trace()
	effarr_mask=effarr[effmask] # apply mask
	timearr_mask=timearr[effmask] # apply mask
	
	#pdb.set_trace()
	
	plt.subplot(1,2,2)
	plt.plot(timearr_mask/3600,effarr_mask,marker='s',linestyle='',color='k',markersize=12,fillstyle='none') # plot simulated points	plt.plot(time_new/3600,eff_smoothed,'--b',linewidth=1.3,dashes=(5,6)) #plot interpolated curve
	
	# interpolate only if there are more than 3 points
	time_new=float('NaN')
	eff_smoothed=float('NaN')
	if len(effarr_mask)>3:
		# Spline interpolation of the efficiency curve
		tck=interpolate.splrep(timearr_mask,effarr_mask,s=0,k=3) # spline interpolation coefficients (order k)
		time_new=np.arange(timearr_mask[0],timearr_mask[-1],300) # create a finer time vector
		effarr_interp=interpolate.splev(time_new,tck,der=0) #evaluate the interpolate curve
	
		# smooth interpolated curve
		wdw=11
		polord=3
		eff_smoothed=signal.savgol_filter(effarr_interp,wdw,polord) #window length wdw, poly order polord
	
		plt.plot(time_new/3600,eff_smoothed,'--b',linewidth=1.3,dashes=(5,6)) #plot interpolated curve
		
		plt.legend(['Simulated\npoints','Guide to\nthe eye'],loc='upper right',prop={'size':9})
	
	plt.rcParams.update({'font.size':16})
	plt.subplots_adjust(left=0.15,top=0.98,wspace=0.4) # adjust space between subplots
	
	plt.ylabel("Efficiency (%)")
	plt.xlabel("Time (h)")
	#plt.ylim(-35, 0)
	#plt.xlim(0, 0.65)
	
	timepts=SimTimePts(results) # Time points that have been run by the simulations
	
	#plt.show()
	return results,timearr_mask,effarr_mask,time_new,eff_smoothed
	
	
def plotcurves():
	[Idark,Vdark]=analyzedat("dark")
	[Ilight,Vlight]=analyzedat("light")
	plt.plot(Vdark,Idark,label='dark')
	plt.hold(True)
	plt.plot(Vlight,Ilight,label='light')
	plt.rcParams.update({'font.size': 20}) # increase fontsize
	plt.legend(loc='upper left')
	plt.show()
	
def findeff(V,I):
# Function returning the efficiency of the cell assuming 1 sun illumination
	
	I=-I # to work with positive currents
	P=V*I
	# Find max powerpoint
	P_L=list(P)
	pm=max(P_L)
	imax=P_L.index(pm)
	
	#[k for k,j in enumerate(P) if j==pm]
	#pdb.set_trace()
	Vmax=V[imax]
	Imax=I[imax]
	
	nu=Imax*Vmax*10/1e3*100 # I is in mA/cm2, corresponding to 10 A/m2
	
	return nu

# function to plot Na profile
def ploth5(h5file="FOR_newNaprofile.h5", folderpath="/home/linux/ieng6/na299x/na299x/DB/Guillaume/Solar_cell_Al_BSF"):
	
	#Open h5 file (file for full stack file)
	hf		= h5py.File(h5file, 'r')
	time 	= hf.get('time')
	ct		= hf.get('si/concentration')
	# ct_10=hf.get('si/concentration/ct_10')
	x		= hf.get('si/x')
	
	pdb.set_trace()
	plt.ion()
	fig,ax=plt.subplots()
	for t in enumerate(time):
		ct_field='ct_'+str(t[0]) # do not use all curve since they are high resolution
		#pdb.set_trace()
		plt.plot(x-x[0],ct[ct_field])
		plt.ylim(1e10,1e16)
		plt.xlim(0,0.8)
		plt.yscale('log')
		
	for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
		ax.get_xticklabels() + ax.get_yticklabels()):
		item.set_fontsize(19)
		
	
	hf.close()
	
	# creates a list of length L filled with NaN
def NaNlist(L): 
	list=np.zeros(L,dtype=np.float) #create array of floats because an array of int cannot be filled with NaN
	list.fill(np.nan) # fill with NaN
	list=list.tolist() #convert to list
	
	return list
	
# Function returning the indices of the time points that have been simulated in the result list output from batchanalysis.py
def SimTimePts(results):
	tme=results["time"]
	tme_array=np.asarray(tme)
	tme_mask=np.isfinite(tme_array)
	timepts=np.where(tme_mask==True)
	
	return timepts[0] # array
	
# Function to make video of efficiency as a function of time
# Takes in argument the images generated at each time
## not currently working, see module videoPID.py
def makePIDvideo(time_new,eff_smoothed):

	#Writer=manimation.FFMpegWriter(fps=30)
	#['ffmpeg']
	#writer=Writer(fps=20, metadata=dict(title='PID simulated efficiency', artist='G.',bitrate=1800)) #record at 20 fps
	#pdb.set_trace()
	# create figure
	figPID,ax=plt.subplots()
	line,=ax.plot(time_new,eff_smoothed)
	plt.xlim(0,12)
	plt.ylim(0,19)
	
	def animate_PID(i,time_new,eff_smoothed,line):
		line.set_data(time_new[:i],eff_smoothed[:i])
		return line,
	
	plt.xlabel('Time (hrs)')
	plt.ylabel('Efficiency (%)')
	
	#plt.show()
	#pdb.set_trace()
	ani = manimation.FuncAnimation(figPID, animate_PID, len(time_new), fargs=[time_new,eff_smoothed,line], blit=False)
	ani.save('Efficiencyplot.mp4', writer='ffmpeg')
	# import images and make video

# Function used to create a log of the parameters and files used in the simulations
def createlog(batchdir,Temp,mseg,clathrate_file, h5file, startstep, endstep, skipNB, sdetemplate, sdevicetemplate):
	# Also add the name of the optical generation file used? (although it is in the sde file)
	arguments=locals() # get all function arguments as a dictionary
	now=datetime.now()
	filename=batchdir+"/AA_logfile_PID_"+now.strftime("%Y%d%m_%H_%M_%S")+".txt"
	fid=open(filename,"w+")
	
	fid.write("Simulation starting time:				"+now.strftime("%Y/%d/%m %H:%M:%S")+"\n")
	fid.write("\n")
	fid.write("Temperature:				"+str(Temp)+" °C\n")
	fid.write("Shunt segregation coefficient:				"+str(mseg)+"\n")
	fid.write("Clathrate conductivity file:				"+clathrate_file+"\n")
	fid.write("Sodium profiles from h5py file:				"+h5file+"\n")
	fid.write("Simulation starting step:				"+str(startstep)+"\n")
	fid.write("Simulation ending step:				"+str(endstep)+"\n")
	fid.write("Step used to skip sodium profiles in the h5py file:				"+str(skipNB)+"\n")
	fid.write("Sentaurus editor template file:				"+sdetemplate+"\n")
	fid.write("Sentaurus device template file:				"+sdevicetemplate+"\n")

# Function replacing expression "expr_search" by "expr_replace" in file "templatename.ext" (in the base directory)
# and saving the updated template "templatename_t<time>.ext" in directory "newfolderpath".
def replace_line(templatename,ext, expr_search,expr_replace,newfolderpath,time,identifier):
# templatename:			name of the template file, eg "sdetemplate"
# extension:			extension of the template file, eg "cmd"
# expr_search:			name of the expression to be replaced, e.g. "nodnum1"
# expr_replace:			name of the new line, e.g. newnodename as newnodename="n_t"+str(int(time))
# newfolderpath:		path of the directory where the simulation data is saved. e.g. "testdir2"
# time:					time of the Na profile (int or float)
# identifier:			type of line being searched. eg: "mesh line"
 
# open sdevice file to find number of lines
	with open(templatename+"."+ext,'r') as fp:
		count_temp = len(fp.readlines(  ))
		fp.close
	
		linelist=[0]*count_temp # preallocate list memory
		newlinelist=[0]*count_temp # preallocate list memory for file with modified line


	with open(templatename+"."+ext,'r') as fp:
		k=0
		successdevice=0
		# save each line into a list
		while 1:
			dataline=fp.readline()
			if dataline=="": # check if the string is empty, meaning end of file
				break
			linelist[k]=dataline # create of list containing one line at each index
			newlinelist[k]=dataline # copy the file into a new list array
			# find line where mesh is defined and change the name of the mesh file
			if re.search(expr_search,linelist[k]): # if the mesh definition is found
				print(identifier + " found")
				newlinelist[k]=re.sub(expr_search,expr_replace,linelist[k]) # replace the nodenumber "nodnum1" in the template by n_t1, n_t2, etc.
				successdevice=successdevice+1
				#pdb.set_trace()
			k=k+1
		fp.close()
		
		if successdevice==0:
			print("*********\n"+identifier + " not found\n**********")
		
		# save the files under a different name (_t0, _t1, _t2, etc)
		# file path
		filepath=os.path.join(newfolderpath,templatename) # path to save file in the right simulation folder

		new_datafile_name=filepath+"_t"+str(int(time))+"."+ext
		f=open(new_datafile_name,"w+")
		for i in range(len(newlinelist)):
			f.write(newlinelist[i])
		f.close()
		print("File "+new_datafile_name+" created.\n")
		
		return new_datafile_name
		
		
