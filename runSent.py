# runSent.py
# Module containing OS interfacing functions

#import subprocess
#from subprocess import Popen, PIPE
import os
import re # Regular Expression package
import pdb # debugger
import numpy as np
import matplotlib.pyplot as plt # plotting package
#import readh5conc
import csv
import h5py

# run svisual to vizualize the defined structure
def run_svisual(meshfile="nodnum1_msh.tdr"):
	# meshfile= "./testdir/
	#run_svisual_os="svisual nodnum1_light_000017_des.tdr" # plot sdevice results
	run_svisual_os="svisual " + meshfile # plot sde created device
	os.system(run_svisual_os)

def run_sde(filename="sde_dvs.cmd"):
# Move to the correct directory where the files are saved
#DB/Guillaume/Solar_cell_Al_BSF/Al_BSF

# Create a character string containing the instructions
#With subprocess the name of the program must be called first
#prog="sde"
#run_sde=" -e -l sde_dvs.cmd"

#It seems that subprocess does not allow to run Sentaurus in batch mode, ie with the option -e (opens GUI instead)
#subprocess.call(["ls","-l"])
#process=Popen([prog,run_sde])

	run_sde_os="sde -e -l " + filename
	
	# Run using the old command
	os.system(run_sde_os) # run Sentaurus Structure Editor
	
def run_sdevice(filename): # filename is the name without extension, ie 'sdevice_dark_des'
	run_sdevice_os="sdevice "+filename
	os.system(run_sdevice_os) # run Sentaurus Device
