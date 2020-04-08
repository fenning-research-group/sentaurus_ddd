** Solar cell device file
** Guillaume von Gastrow

File
{
	**** INPUT FILES
	* geometry, contacts, doping and mesh
	Grid ="nodnum1_msh.tdr"
	* physical parameters
	Parameter = "sdevice.par"
	* IlluminationSpectrum= "@pwd@/spectrum/am15g_1.2um.txt"
	* Import the optical generation input in the mesh file, created in the structure editor file (section optical generation profile)
	OpticalGenerationInput= "nodnum1_msh.tdr"
	**** OUTPUT FILES
	* to visualize distributed variables
	Plot = "nodnum1_light_des.tdr"
	* to visualize electrical characteristics at the electrodes
	Current= "nodnum1_light_des.plt"

}

Electrode
{
	* defines which contacts have to be treated as electrodes
	* & initial boundary conditions
	* obviously, electrode names must match the contact names of the
	* sde_dvs.cmd file
	{ name="em_contact" voltage=0.0 }
	{ name="base_contact" voltage=0.0 }
}

# Start Physics section
Physics
{
	Mobility (
		DopingDependence
			)
	Recombination (
					SRH
					)

	* Define area of the contacts (in um)
	AreaFactor=2e7

	* Use the optical generation file imported and added to the mesh in the Sentaurus structure Editor file
	*Optics(
	*    OpticalGeneration(
	*      ReadFromFile(
	*       TimeDependence(
	*         WaveTime= (1, 2)
	*         WaveTSlope= 0.05
	*        )
	*       )
	*    )
	*  )
	* }
  
	Optics (
		OpticalGeneration (
			SetConstant (
			Value = 2e21
			)
		)
	)
} # End of main Physics section

* Include surface recombination at the silicon/SiNx interface
* SRH recombination parameters will be defined separately below
# Physics (RegionInterface="SiNx-region/Si-profile-region") { Recombination(surfaceSRH) }


Plot{
OpticalGeneration
}

CurrentPlot{
OpticalGeneration(Integrate(Semiconductor) )
  OpticalGeneration(Integrate(material="Silicon") )
  }

Plot
{
* On-mesh-defined variables to be saved in the .tdr output file
*- Doping Profiles
		Doping DonorConcentration AcceptorConcentration
*- Charge, field, potential and potential energy
		SpaceCharge
		ElectricField/Vector Potential
		BandGap EffectiveBandGap BandGapNarrowing ElectronAffinity
		ConductionBandEnergy ValenceBandEnergy
*- Carrier Densities:
		EffectiveIntrinsicDensity IntrinsicDensity
		eDensity hDensity
		eQuasiFermiEnergy hQuasiFermiEnergy
*- Currents and current components:
		eGradQuasiFermi/Vector hGradQuasiFermi/Vector
		eMobility hMobility eVelocity hVelocity
		Current/Vector eCurrent/Vector hCurrent/Vector
		eDriftVelocity/Vector hDriftVelocity/Vector
*- SRH & interfacial traps
SRHrecombination
tSRHrecombination
*- Band2Band Tunneling & II
		eBand2BandGeneration hBand2BandGeneration Band2BandGeneration
		eAvalanche hAvalanche Avalanche
		
}
Math
{
* use previous two solutions (if any) to extrapolate next
Extrapolate
* use full derivatives in Newton method
Derivatives
* control on relative and absolute errors
-RelErrControl
* relative error= 10^(-Digits)
Digits=5
* absolute error
Error(electron)=1e8
Error(hole)=1e8
* numerical parameter for space-charge regions
eDrForceRefDens=1e10
hDrForceRefDens=1e10
* maximum number of iteration at each step
Iterations=20
* choosing the solver of the linear system
Method=ParDiSo

* display simulation time in 'human' units
Wallclock
* display max.error information
CNormPrint
* to avoid convergence problem when simulating defect-assisted tunneling
NoSRHperPotential
}

* Solve the transport equations using the quasistationary command
Solve
{
	* EQUILIBRIUM 
	*coupled {poisson} * do not plot as we are not interested in equilibrium conditions
	
	** TURN-ON
	** if needed, add a transient step to help solving due to the large amount of optical generated carriers
	Transient (
    InitialTime= 0 FinalTime= 1.2
    InitialStep= 1 MaxStep= 1 MinStep= 1e-5
  ) { Coupled {Poisson Electron Hole} }
 
	* decreasing em_contact to goal
	* * quasistationary (InitialStep = 0.010 MaxStep = 0.050 MinStep=0.001
					* * Goal {name= "em_contact" voltage = -0.1}
					* * plot { range=(0, 1) intervals=1 }
					* * )
					* * {coupled {poisson electron hole} }

	** raising em_contact to goal
	** negative part
	quasistationary (InitialStep = 0.01 MaxStep = 0.050 MinStep=0.005 * Max voltage step of 40 mV in the IV curve
	Goal {name= "em_contact" voltage = 0}
	)
	{coupled {poisson electron hole} }
	
	quasistationary (InitialStep = 0.005 MaxStep = 0.050 MinStep=0.001
	Goal {name= "em_contact" voltage = 0.8}
	plot { range=(0, 1) intervals=2 }
	)
	{coupled {poisson electron hole} }
}