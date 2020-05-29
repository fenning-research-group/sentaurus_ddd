import os
from string import Template
from datetime import datetime
import subprocess
from typing import Union, List
import pid.utils as utils
import logging
import numpy as np
import h5py
import pandas as pd
import shutil
from scipy import interpolate

na_profile_type = np.dtype([('x (um)', 'd'), ('c (1/cm3)', 'd')])
jv_dtype = np.dtype([('voltage (V)', 'd'), ('current (mA/cm2)', 'd')])
na_file = "two_layers_D1=4E-16cm2ps_D2=1E-14cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-10_m1.0e+00_pnp.h5"


class PIDModel:
    """
    This class provides methods to run the electrical module of the DDD model
    """
    _NaProfileFileName: str = None
    _batch_band_diagram_template: str = 'batch_band_diagram.tcl'
    _cSi: float = 5e22
    _dshuntname: str = "dshunt1"
    _illuminated: bool = True
    _logger: logging.Logger = None
    _loggingLevels: dict = {'NOTSET': logging.NOTSET,
                            'DEBUG': logging.DEBUG,
                            'INFO': logging.INFO,
                            'WARNING': logging.WARNING,
                            'ERROR': logging.ERROR,
                            'CRITICAL': logging.CRITICAL}
    _sde_template_3d: str = 'metallic_shunt_3d_sde_dvs.scm'
    _sde_template_2d: str = 'metallic_shunt_sde_dvs.scm'
    _sde_template_traps: str = 'shunt_traps_sde_dvs.scm'
    _sde_template_traps_cylindrical: str = 'shunt_traps_sde_dvs_cylindrical.scm'
    _sde_template_3_cylindrical: str = 'metallic_shunt_cylindrical_sde_dvs.scm'
    _sdevice_template_2d: str = 'sdevice_light_des.cmd'
    _sdevice_template_3d: str = 'sdevice_light_3d_des.cmd'
    _simulation_time: np.ndarray = np.array([])
    _simulation_indices: np.ndarray = np.array([])
    _spectrum: str = 'am15-g_reduced.txt'
    _parfile_template: str = 'sdevice.txt'

    def __init__(self, folder_path: str = "/home/linux/ieng6/na299x/na299x/DB/Erick/PID/test",
                 na_h5_file: str = na_file,
                 segregation_coefficient: float = 10,
                 use_srv: bool = False,
                 clathrate_file: str = "clathrates_cond.csv",
                 contact_length: float = 5,
                 cell_width: float = 50,
                 cell_length: float = 50,
                 optical_generation: str = None):
        """

        Parameters
        ----------
        folder_path: str
            Name of the directory where the generated sde and sdevice files will be saved. Example: "./test_dir"
        na_h5_file: str
            Name of the h5py file containing the sodium migration profiles at each time point. eg "FOR_newNaprofile.h5"
        segregation_coefficient: [int, float]
            Value of the segregation coefficient of sodium from Si bulk to the stacking fault. (Default 10)
        use_srv: bool
            True if the model incorporates SRV, false otherwise
        clathrate_file: str
            Name of the file containing the fit of clathrate resistivity as a function of Na to Si ratio
        contact_length: float
            The length of the metal contact
        cell_width: float
            The total length of the cell
        optical_generation: str
            The path to the optical generation file
        """

        self._clathrateFile: str = clathrate_file
        self._cwd: str = os.path.dirname(__file__)
        self._folder_name: str = folder_path
        self._folder_path: str = folder_path
        self._NaProfileFileName: str = na_h5_file
        self._segregation_coefficient: float = segregation_coefficient
        self._useSRV: bool = use_srv
        self._contact_length = contact_length
        self._cell_width = cell_width
        self._cell_length = cell_length
        if optical_generation is None:
            self._optical_generation_file = os.path.join(self._cwd, 'input_optical_generation.plx')
        else:
            self._optical_generation_file = optical_generation

        with h5py.File(na_h5_file, 'r') as hf:
            group_time = hf['time']
            self._na_profile_time_array = np.array(hf['time'])
            self._temperature = float(group_time.attrs['temp_c'])

    def run_pid(self, start_step: int = 0, end_step: int = 0, skip_nb: int = 0, overwrite_folder: bool = False):
        """
        Start the sdevice simulation for each C(t) in [startstep:skip_nb:endstep]

        Parameters
        ----------
        start_step: int
            Number of the step where simulations will start. Should be 0 unless a different starting step is wanted
            (for instance if previous steps have already been run before).
        end_step: int
            Number of the step where the simulations will end. If set to 0, simulations will run at all time points in
            the h5file containing the sodium migration profiles
        skip_nb: int
            Number of time points to skip in the sodium migration dataset at each iteration (0 by default).
            If set to 1, every second point will be run. If set to 2, every third point. etc.
        overwrite_folder: bool
            If the directory already exists, False will create a new one while True will save in the same one.
            By default changedir=False.

        Returns
        -------

        """
        skip_nb = abs(int(skip_nb))
        start_step = abs(int(start_step))
        end_step = abs(int(end_step))

        if skip_nb == 0:
            self.log('Trying to set skip_nb={0:d} in run_pid method. Reverting to skip_nb=1.'.format(skip_nb),
                     level='WARNING')
            skip_nb = 1
        elif start_step > end_step:
            self.log('Trying to set start_step={0:d} > end_step{1:d} in run_pid method. Reverting to skip_nb=1.'.format(
                start_step, end_step),
                level='WARNING')
            start_step = 0
            end_step = 0
            skip_nb = 1

        self._folder_path: str = self.prepare_folder(self._folder_name, overwrite_folder)
        self._logger: logging.Logger = self.__create_logger(logger_id='pid_logger')

        # Start simulation log
        now = datetime.now()
        # time_stamp = now.strftime('%Y%m%d')
        iso_date = now.astimezone().isoformat()
        self.log('Simulation starting {0}'.format(iso_date))
        self.log('Temperature: {0} °C'.format(self._temperature))
        self.log('Shunt segregation coeffiecient: {0:g}'.format(self._segregation_coefficient))
        self.log('Sodium profiles simulation file: {0}'.format(self._NaProfileFileName))
        self.log('Simulation starting step: {0}'.format(start_step))
        self.log('Step used to skip sodium profiles in the h5py file: {0}'.format(skip_nb))
        self.log('Simulation ending step: {0}'.format(end_step))

        # factor=50 # factor = 100 to have a final depth of 600 nm, as in "For_JV_85C.h5" the final depth is 6 nm.
        # NOTE: Later, with the correct Na profiles, it should be just one.
        factor = 1
        with h5py.File(self._NaProfileFileName, 'r') as hf:
            time = np.array(hf.get('time'))
            depth = np.array(hf.get('L2/x'))

        x = (depth - depth[0]) * factor  # Difference because x does not start at 0
        self._simulation_time = np.array([time[i] for i in range(start_step, end_step + skip_nb, skip_nb)])
        self._simulation_indices = np.array([i for i in range(start_step, end_step + skip_nb, skip_nb)])

        for i in range(start_step, end_step + skip_nb, skip_nb):
            # Get the Na concentration at time t = ti
            ct_ds = 'ct_{0:d}'.format(i)
            with h5py.File(self._NaProfileFileName, 'r') as hf:
                group_si_c = hf['/L2/concentration']
                ct = np.array(group_si_c[ct_ds])
            # Estimate the conductivity
            c_shunt = ct * self._segregation_coefficient
            # Save conductivity profile in a .plx file
            concentration_filename = "na_profile_t{0:.0f}.plx".format(time[i])
            concentration_fullname = os.path.join(self._folder_path, concentration_filename)
            with open(concentration_fullname, 'w') as fp:
                fp.write("\"DeepLevels\"\n")
                line = "75E-3\t0\n"
                line += "50E-3\t0\n"
                line += "25E-3\t0\n"
                line += "1E-4\t0\n"
                fp.write(line)
                for xi, c in zip(x, c_shunt):
                    line = "{0:1.4E}\t{1:1.4E}\n".format(xi, c)
                    fp.write(line)
            self.log('Created concentration file \'{0}\'.'.format(concentration_fullname))
            # Update the parfile with the concentration
            self.srv_param(cNa=ct, time=time[i], use_srv=True)
            [sde_filename, sdevice_filename] = self.update_files(concentration=ct, x=x, time=time[i])
            # Run sde
            success: bool = self.run_sde(input_file=sde_filename)
            if not success:
                raise RuntimeError('Error running sde. Stopping program execution.')
            # Run sdevice
            self.run_sdevice(input_file=sdevice_filename)

        self.log('Converting JV plots to tdr files...')
        self.plt2tdr()
        final_string = "Results in {0}".format(self._folder_path)
        self.log(final_string)

    def request_time_points(self, t_min: float, t_max: float, t_steps: int, spacing: str = "log"):
        """
        Define the simulated points to run. This function generates a list of indices corresponding to times at which
        the simulated Na profiles are drawn from the h5 file.

        Parameters
        ----------
        t_min: float
            The minimum time of the simulation (s)
        t_max: float
            The maximum time of the simulation (s)
        t_steps: int
            The number of time steps to simulate
        spacing: List[str]
            The type of spacing (linear, log, geometric)

        Returns
        -------
        np.ndarray
            A list with the indices corresponding to the desired samplig
        """
        min_delta = self._na_profile_time_array[1] - self._na_profile_time_array[0]

        if spacing not in ["linear", "log", "geometric"]:
            spacing = "log"  # Default to log

        if spacing == "linear":
            proposed_times = np.linspace(t_min, t_max, t_steps + 1)
        elif spacing == "log":
            t_min = t_min if t_min > 0 else 1E-1
            proposed_times = np.logspace(np.log10(min_delta), np.log10(t_max), t_steps)
            proposed_times = np.insert(proposed_times, 0, 0.0)
        else:
            proposed_times = utils.geometric_series_spaced(
                max_val=t_max,
                min_delta=min_delta,
                steps=t_steps
            )

        requested_indices = utils.get_indices_at_values(
            x=self._na_profile_time_array,
            requested_values=proposed_times
        )
        return requested_indices

    def run_pid_metal(self, requested_indices: Union[List[int], np.ndarray], overwrite_folder: bool = False,
                      cylindrical: bool = False, activated_na: float = 1.0, fixed_conductivity: float = 0,
                      three_dimensional: bool = True,
                      shunt_depth: float = 0, metal_work_function: float = 4.05,
                      uniform_profile: bool = False,
                      use_srv: bool = False, illuminated: bool = True):
        """
        Start the sdevice simulation for each C(t) in [startstep:skip_nb:endstep]

        Parameters
        ----------
        requested_indices: Union[List[int], np.ndarray]
            An array with the indices corresponding to Na profiles as a function of time
        overwrite_folder: bool
            If the directory already exists, False will create a new one while True will save in the same one.
            By default changedir=False.
        cylindrical: bool
            If true, apply cylindrical symmetry to the structure.
        three_dimensional: bool
            True if running a 3D simulation
        activated_na: float
            The fraction of sodium that is ionized [0 - 1]
        fixed_conductivity: float
            If > 0 used a fixed value of the conductivity for the whole shunt (Ohm cm)^{-1}
        shunt_depth: float
            If 0, ignore. If > 0 use it as constant shunt depth (um)
        metal_work_function: float
            The work function of the metal shunt (default 4.05 eV)
        uniform_profile: bool
            True if using a uniform conductivity in the shunt based on the concentration
        use_srv: bool
            True if estimating SRV
        illuminated: bool
            If set to True, simulates an illuminated device, otherwise simulates a dark I-V. Default: True

        Returns
        -------

        """
        self._illuminated = illuminated
        self._folder_path: str = self.prepare_folder(self._folder_name, overwrite_folder)
        self._logger: logging.Logger = self.__create_logger(logger_id='pid_logger')

        # Start simulation log
        now = datetime.now()
        # time_stamp = now.strftime('%Y%m%d')
        iso_date = now.astimezone().isoformat()
        self.log('Simulation starting {0}'.format(iso_date))
        self.log('Temperature: {0} °C'.format(self._temperature))
        self.log('Shunt segregation coeffiecient: {0:g}'.format(self._segregation_coefficient))
        self.log('Sodium profiles simulation file: {0}'.format(self._NaProfileFileName))

        with h5py.File(self._NaProfileFileName, 'r') as hf:
            time = np.array(hf.get('time'))
            depth = np.array(hf.get('L2/x'))

        x = (depth - depth[0])
        time_msk = np.zeros(len(time), dtype=bool)
        time_msk[requested_indices.astype(int)] = True
        self._simulation_time = time[time_msk]
        self._simulation_indices = requested_indices
        requested_times = time[time_msk] / 3600
        requested_times_str = ['{0:.3f}'.format(rt) for rt in requested_times]
        self.log('Simulation time points : {0} hr'.format(', '.join(requested_times_str)))

        shunt_depth = shunt_depth if shunt_depth > 0 else np.amax(x)

        for i in requested_indices.astype(int):
            # Get the Na concentration at time t = ti
            ct_ds = 'ct_{0:d}'.format(i)
            with h5py.File(self._NaProfileFileName, 'r') as hf:
                group_si_c = hf['/L2/concentration']
                ct = np.array(group_si_c[ct_ds])
            # Estimate the conductivity
            conductivity = self.conductivity_model(ct[:], activated_na=activated_na)
            if uniform_profile:
                idx = (np.abs(conductivity - conductivity[0] / 100)).argmin()
                conductivity = conductivity[idx]
                [sde_filename, sdevice_filename] = self.input_files_metal(concentration=ct, x=x, time=time[i],
                                                                          three_dimensional=False,
                                                                          cylindrical=cylindrical,
                                                                          activated_na=activated_na,
                                                                          fixed_conductivity=conductivity,
                                                                          shunt_depth=shunt_depth,
                                                                          metal_work_function=metal_work_function,
                                                                          use_srv=use_srv)
                self.log('Conductivity = {0:.3E} S/cm, x = {1:.3E} um, [Na] = {2:.3E} cm^-3'.format(
                    conductivity, x[idx], ct[idx]
                ))
            else:
                # # Save conductivity profile in a .plx file
                # conductivity_filename
                # = "conductivity_t{0:.0f}.plx".format(time[i])
                # conductivity_fullname = os.path.join(self._folder_path, conductivity_filename)
                # with open(conductivity_fullname, 'w') as fp:
                #     fp.write("\"MetalConductivity\"\n")
                #     line = "75E-3\t0\n"
                #     line += "50E-3\t0\n"
                #     line += "25E-3\t0\n"
                #     line += "1E-4\t0\n"
                #     fp.write(line)
                #     for xi, c in zip(x, conductivity):
                #         line = "{0:1.4E}\t{1:1.4E}\n".format(xi, c)
                #         fp.write(line)
                # self.log('Created concentration file \'{0}\'.'.format(conductivity_fullname))
                # Update the parfile with the concentration
                # self.srv_param(cNa=ct, time=time[i], use_srv=True)
                self.log('Processing time t = {0:.0f} s'.format(time[i]))

                [sde_filename, sdevice_filename] = self.input_files_metal(concentration=ct, x=x, time=time[i],
                                                                          three_dimensional=three_dimensional,
                                                                          cylindrical=cylindrical,
                                                                          activated_na=activated_na,
                                                                          fixed_conductivity=fixed_conductivity,
                                                                          shunt_depth=shunt_depth,
                                                                          metal_work_function=metal_work_function,
                                                                          use_srv=use_srv)
            # Run sde
            success: bool = self.run_sde(input_file=sde_filename)
            if not success:
                raise RuntimeError('Error running sde. Stopping program execution.')
            # Run sdevice
            self.run_sdevice(input_file=sdevice_filename)

        self.log('Converting JV plots to tdr files...')
        self.plt2tdr()
        final_string = "Results in {0}".format(self._folder_path)
        self.log(final_string)

    def conductivity_model(self, concentration: np.ndarray, activated_na: float = 1) -> np.ndarray:
        """
        Implementation of the conductivity_model model.

        =====================
        Model simplifications
        =====================
        1. The Na to Si ratio in the stacking fault is obtained from the ratio between Na concentration and Si
            concentration in the bulk of a perfect crystal (does not consider the specific geometry of a stacking fault)
        2. Conductivity is calculated based on depth-resolved Hall-effect measurements of mobility and carrier density
            in Na-implanted Si (Korol et al)
           Reference:
            Korol, V. M. "Sodium ion implantation into silicon." Physica status solidi (a) 110.1 (1988): 9-34.

        Parameters
        ----------
        concentration: np.ndarray
            The sodium concentration in the Si bulk
        activated_na: float
            The fraction of Na that is activated

        Returns
        -------
        np.ndarray
            The conductivity_model profile
        """

        # Na concentration in the shunt
        cshunt = concentration * self._segregation_coefficient * activated_na

        # Clathrate model, not used because not realistic at our Na concentrations
        # TODO: consider removing the clathrate model
        # if False:  # Skip this section with a trivial false condition
        #     # Model for Na density in shunt: ratio between Na density in shunt and Si bulk density in a perfect
        #     crystal
        #
        #     Na_Si_ratio = cshunt / cSi
        #     # Add condition [Na]/[Si]=1 if higher than Si density
        #
        #     # Later import the coefficients directly from the conductivity_model .csv file
        #     # coefficients at 60 deg C
        #     # a=24.960488582173806
        #     # b=-1.7580985174592794
        #
        #     # coefs 70 C
        #     a = 24.496024405605336
        #     b = -1.681009278808291
        #
        #     # Calculate conductivity_model profile
        #     sigma = 10 ** (a * Na_Si_ratio + b)  # S/cm

        # Model based on implantation data
        # Korol, V. M. "Sodium ion implantation into silicon." Physica status solidi (a) 110.1 (1988): 9-34.
        # Fitting of coefficients in Extract_NaImp.py
        coord = -11.144769029961262
        slope = 0.717839509854622

        sigma = (10 ** coord) * (cshunt ** slope)  # S/cm

        return sigma

    def update_files(self, concentration: np.ndarray, x: np.ndarray,
                     time: Union[int, float]) -> List[str]:
        """
        This method modifies the Sentaurus input files to include updated shunt depth and modified external files
        (including the .par file). Creates a shunt with spatially varying conductivity_model based on doping profiles
        calculated at each time.

        Parameters
        ----------
        concentration: np.ndarray
            Sodium concentration from the h5 file
        x: np.ndarray
            depth from the h5 file (list)
        time: Union[int, float]
            time from the h5 file (int or float)
        Returns
        -------
        List[str]
            [sde_file_path, sdevice_file_path]
        """
        # The directory to save the data to
        output_folder = self._folder_path
        # name for the mesh file
        mesh_name = "n_t{0:d}".format(int(time))
        # Calculate segregation coefficient at each depth and the depth of the shunt
        cseg = np.abs(concentration) * self._segregation_coefficient
        # pdb.set_trace()
        # k = 0
        # # ATTENTION: need to modify this as the Na concentration is not a meaningful value here. The shunt depth is
        # # highly sensitive to the conductivity_model model and to segregation_coefficient.
        # while cseg[k] > 0.1 and k < L:
        #     # Consider only the part of the [Na] profile in the stacking fault that is at least 1 cm-3
        #     # (to remove useless extra points). TO CHANGE.
        #     k += 1
        idx_threshold = cseg >= 0.1
        # cseg = cseg[idx_threshold]
        x = x[idx_threshold]

        # **************************************************************************************************************
        # ATTENTION! This should be much more reliable than using cseg[k]>0.1, to be tested.
        # Or maybe we shouldn't limit the shunt depth at all.
        # while(PIDmodel.conmodel(cseg[k],60,mseg) > 1e-3 and k<L-1): # Limit the shunt depth to conductivities higher
        # than 1e-3 S/cm (assuming 1e-3 is lower than the limit conductivity_model for shunting)
        # k=k+1
        # **************************************************************************************************************
        if time == 0 or len(x) == 0:
            shunt_depth = 0
        elif len(x) == 1:
            shunt_depth = x[0]
        else:
            shunt_depth = np.amax(x)  # Depth of the profile in um. The depth does not start at 0 so subtract x[0]

        # Define an arbitrary shunt depth in case of very low concentrations
        # that would lead to a shunt depth of 0 um and cause Sentaurus to crash
        if abs(shunt_depth) < 1E-4:
            shunt_depth = 0.1

        # NOTE: CHECK IF IT MATTERS IF THE EXTERNAL PROFILE IS DEEPER THAN THE DEPTH DEFINED HERE. CURRENTLY IT IS.
        # *************** Generate the sde file based on the template *******************
        # Prepare the substitutions
        conductivity_filename = "na_profile_t{0:.0f}.plx".format(time)
        conductivity_fullname = os.path.join(self._folder_path, conductivity_filename)
        mesh_node_file = os.path.join(self._folder_path, mesh_name)
        substitutions_sde = {
            'dshunt_name': self._dshuntname,
            'shunt_depth': shunt_depth,
            'shunt_profile_file': conductivity_fullname,
            'optical_generation_file': self._optical_generation_file,
            'nodenum': mesh_node_file,
            'contact_length': self._contact_length,
            'device_length': self._cell_width,
            'c0': np.amax(concentration * self._segregation_coefficient)
        }
        # Load the template file
        template_file_sde = open(os.path.join(self._cwd, self._sde_template_2d), 'r')
        src = Template(template_file_sde.read())
        # perform the substitutions
        result = src.substitute(substitutions_sde)
        template_file_sde.close()
        output_filename_sde = 'sde_dvs_t{0:d}.cmd'.format(int(time))
        output_filename_sde = os.path.join(self._folder_path, output_filename_sde)
        # Save the sde file
        with open(output_filename_sde, 'w') as output_file:
            output_file.write(result)

        self.log('Created Sentaurus input file \'{0}\'.'.format(output_filename_sde))

        # *************** Generate the sdevice file based on the template *******************
        # Prepare the substitutions
        node_name_msh = '{0}_msh'.format(os.path.join(output_folder, mesh_name))
        if self._illuminated:
            out_mode = 'light'
        else:
            out_mode = 'dark'
        node_name_out = '{0}_{1}_des'.format(os.path.join(output_folder, mesh_name), out_mode)
        parfile = '{0}_t{1}.par'.format(os.path.join(output_folder, 'sdevice'), int(time))
        substitutions_sdevice = {
            'nodenum': node_name_msh,
            'parfile': parfile,
            'node_out': node_name_out,
            'illumination_spectrum': os.path.join(self._cwd, self._spectrum),
            'area_factor': '{0:.4e}'.format(1E11 / (25 * self._cell_width)),  # 2D
            'folder_path': self._folder_path,
            'dshunt_name': self._dshuntname
        }
        # 'area_factor': '{0:.4e}'.format(1E11 / (self._device_length)), # 2D
        # Load the template file
        template_file_sdevice = open(os.path.join(self._cwd, self._sdevice_template_2d), 'r')
        src = Template(template_file_sdevice.read())
        # perform the substitutions
        result = src.substitute(substitutions_sdevice)
        template_file_sdevice.close()
        output_filename_sdevice = 'sdevice_des_t{0}.cmd'.format(int(time))
        output_filename_sdevice = os.path.join(self._folder_path, output_filename_sdevice)
        # Save the sde file
        with open(output_filename_sdevice, 'w') as output_file:
            output_file.write(result)

        self.log('Created Sentaurus input file \'{0}\'.'.format(output_filename_sdevice))

        return [output_filename_sde, output_filename_sdevice]

    def input_files_metal(self, concentration: np.ndarray, x: np.ndarray,
                          time: Union[int, float], use_srv: bool = False,
                          three_dimensional: bool = True,
                          cylindrical: float = False,
                          activated_na: float = 1.0,
                          fixed_conductivity: float = 0,
                          shunt_depth: float = 0,
                          metal_work_function: float = 4.05,
                          n_regions: int = 100) -> List[str]:
        """
        This method modifies the Sentaurus input files to include updated shunt depth and modified external files
        (including the .par file). Creates a shunt with spatially varying conductivity_model based on doping profiles
        calculated at each time.

        Parameters
        ----------
        concentration: np.ndarray
            Sodium concentration from the h5 file
        x: np.ndarray
            depth from the h5 file (list)
        time: Union[int, float]
            time from the h5 file (int or float)
        use_srv: bool
            If true updates the parameter file with the estimated SRV
        three_dimensional: float
            True if a 3D simulation otherwise 2D
        cylindrical: bool
            If true, activate cylindrical symmetry
        activated_na: float
            The fraction of sodium that is ionized.
        fixed_conductivity: float
            If the value is 0, ignore it. If > 0 use the value as a constant conductivity in the shunt (Ohm cm)^{-1}
        shunt_depth: float
            If the value is 0, ignore. If >0 use it as shunt depth (in um)
        metal_work_function: float
            The work function of the metal shunt (default 4.05 eV)
        n_regions: int
            The number of regions to interpolate the conductivity profile into

        Returns
        -------
        List[str]
            [sde_file_path, sdevice_file_path]
        """
        # The directory to save the data to
        output_folder = self._folder_path
        # name for the mesh file
        mesh_name = "n_t{0:d}".format(int(time))
        # Calculate segregation coefficient at each depth and the depth of the shunt
        cseg = np.abs(concentration) * self._segregation_coefficient

        # Get the conductivity
        if fixed_conductivity == 0:
            conductivity = self.conductivity_model(cseg, activated_na=activated_na)
            # idx_threshold = conductivity >= 1E-3
            # conductivity = conductivity[idx_threshold]
            # x = x[idx_threshold]
            if time == 0:
                shunt_depth = 0
                dy = 0
                y = np.zeros(1)
                sigma = np.zeros(1)
                n_regions = 0
                self.log('1. Shunt depth = 0, Variable conductivity, Time = 0, sigma = [0]')
            elif len(x) == 1:
                shunt_depth = x[0]
                dy = shunt_depth
                y = np.zeros(1)
                sigma = conductivity[0]
                self.log(
                    '2. Shunt depth = {0:.3E} um, Variable conductivity,  time = {1:.3f} h, sigma = {2:.3E}'.format(
                        shunt_depth, time / 3600, sigma[0]
                    )
                )
            else:
                shunt_depth = np.amax(x)  # Depth of the profile in um. The depth does not start at 0 so subtract x[0]

                if len(x) >= n_regions:
                    # Interpolate the profile to get only n_regions points
                    ymax = np.amax(x)
                    if 0 < shunt_depth < ymax:
                        ymax = shunt_depth
                    y = np.linspace(0., ymax, n_regions + 1)
                    dy = y[1] - y[0]
                    f = interpolate.interp1d(x, conductivity, kind='slinear', fill_value="extrapolate")
                    y = y[0:-1]
                    sigma = f(y)
                    debug_str = '3. len(x) = {0:d} Shunt depth = {1:.3E} um, Variable conductivity,  ' \
                                'time = {2:.3f} h, sigma = {3}'
                    self.log(
                        debug_str.format(len(x), shunt_depth, time / 3600, np.array_str(sigma))
                    )
                elif len(x) >= 1:
                    if shunt_depth > 0:
                        idx = int((np.abs(shunt_depth - x)).argmin())
                    y = x[0:idx]
                    dy = y[1] - y[0]
                    sigma = conductivity[0:idx]
                    debug_str = '4. len(x) = {0:d} Shunt depth = {1:.3E} um, Variable conductivity,  ' \
                                'time = {2:.3f}, sigma = {3}'
                    self.log(
                        debug_str.format(len(x), shunt_depth, time / 3600, np.array_str(sigma))
                    )
                else:
                    dy = 0
                    y = np.zeros(1)
                    sigma = np.zeros(1)
                    self.log('5. len(x) = 0. Shunt depth = 0, Variable conductivity, Time = 0, sigma = [0]')
        else:
            if time == 0:
                shunt_depth = 0
                dy = 0
                y = np.zeros(1)
                sigma = np.array([0])
                n_regions = 0
                self.log('6. Shunt depth = 0, Fixed conductivity, Time = 0, sigma = [0]')
            else:
                shunt_depth = abs(shunt_depth)
                fixed_conductivity = abs(fixed_conductivity)
                dy = shunt_depth
                sigma = np.array([fixed_conductivity])
                y = np.zeros(1)
                self.log(
                    '7. Shunt depth = 0, Fixed conductivity, Time = {0:.3f}, sigma = {1:.3E}'.format(
                        time / 3600,
                        fixed_conductivity
                    )
                )

        box_indices = ' '.join(['%d' % i for i in range(len(sigma))])

        # *************** Generate the sde file based on the template *******************
        # Prepare the substitutions
        # conductivity_filename = "na_profile_t{0:.0f}.plx".format(time)
        # conductivity_fullname = os.path.join(self._folder_path, conductivity_filename)
        mesh_node_file = os.path.join(self._folder_path, mesh_name)

        # Load the template file
        symmetry: str = '* No symmetry'
        if three_dimensional:
            sde_shunt_regions = ""
            sde_shunt_ref = ""
            if shunt_depth > 0:
                shunt_regions = utils.get_triangle_regions(
                    shunt_depth=shunt_depth, n_regions=n_regions,
                    x_center=(self._cell_width / 2.), y_center=(self._cell_length / 2.)
                )
                sde_region_template = "(sdegeo:create-polygon " \
                                      "(list ${positions}) " \
                                      "\"Metal\" " \
                                      "\"${region_name}\")\n"
                sde_region_template += "(sdegeo:extrude (list (car (find-face-id ${center} ))) 0.005)"
                src = Template(sde_region_template)
                # The base conductivity of
                for reg in shunt_regions:
                    points = reg['points']
                    point_list = ""
                    for i, p in enumerate(points):
                        p_str = "(position {0:.6f} {1:.6f} {2:.6f})".format(p['x'], p['y'], p['z'])
                        point_list += p_str
                        if i < (len(points) - 1):
                            point_list += " "
                    pos_center = reg['center']
                    reg_center_str = '(position {0:.8f} {1:.8f} {2:.8f})'.format(
                        pos_center[0], pos_center[1], pos_center[2]
                    )
                    substitutions = {
                        'region_name': 'shunt.Region.{0:d}'.format(reg['region']), 'positions': point_list,
                        'center': reg_center_str
                    }
                    sde_shunt_regions += src.safe_substitute(substitutions)
                    sde_shunt_regions += "\n"

                sde_shunt_ref = '(sdedr:define-refeval-window ' \
                                '"RefWin.Shunt" "Cuboid" ' \
                                '(position {x0:.6f} {y0:.6f} 0) ' \
                                '(position {x1:.6f} {y1:.6f} (- {z:.6f})) )'
                sde_shunt_ref += "\n"
                sde_shunt_ref += '(sdedr:define-refinement-size "RefDef.Shunt" ' \
                                 '(/ SHUNT_DEPTH 10) (/ SHUNT_DEPTH 10) (/ SHUNT_DEPTH 10) ' \
                                 '(/ SHUNT_DEPTH 10) (/ SHUNT_DEPTH 10) (/ SHUNT_DEPTH 10))'
                sde_shunt_ref += "\n"
                sde_shunt_ref += '(sdedr:define-refinement-placement "PlaceRF.Shunt" "RefDef.Shunt" "RefWin.Shunt")'
                x_center = self._cell_width / 2.
                y_center = self._cell_length / 2.
                sde_shunt_ref = sde_shunt_ref.format(
                    x0=(x_center - shunt_depth), y0=(y_center - shunt_depth), x1=(x_center + shunt_depth),
                    y1=(y_center + shunt_depth), z=shunt_depth
                )

            substitutions_sde = {
                'shunt_depth': shunt_depth,
                'optical_generation_file': self._optical_generation_file,
                'nodenum': mesh_node_file,
                'contact_length': self._contact_length,
                'cell_length': self._cell_length,
                'cell_width': self._cell_width,
                'shunts': sde_shunt_regions,
                'shunt_refinement': sde_shunt_ref
            }
            template_file_sde = open(os.path.join(self._cwd, self._sde_template_3d), 'r')
            area_factor = 1E11 / (self._cell_length * self._cell_width)
        else:
            substitutions_sde = {
                'dshunt_name': self._dshuntname,
                'shunt_depth': shunt_depth,
                'optical_generation_file': self._optical_generation_file,
                'nodenum': mesh_node_file,
                'contact_length': self._contact_length,
                'device_length': self._cell_width,
                'shunt_positions': ' '.join(map(str, y)),
                'shunt_sigmas': ' '.join(map(str, sigma)),
                'shunt_index': box_indices,
                'shunt_dy': dy
            }
            if cylindrical:
                template_file_sde = open(os.path.join(self._cwd, self._sde_template_3_cylindrical), 'r')
                area_factor = 1E11 / (np.pi * (self._cell_width / 2) ** 2)
                symmetry: str = 'Cylindrical'
            else:
                template_file_sde = open(os.path.join(self._cwd, self._sde_template_2d), 'r')
                area_factor = 1E11 / self._cell_width

        src = Template(template_file_sde.read())
        # perform the substitutions
        result = src.substitute(substitutions_sde)
        template_file_sde.close()
        output_filename_sde = 'sde_dvs_t{0:d}.cmd'.format(int(time))
        output_filename_sde = os.path.join(self._folder_path, output_filename_sde)
        # Save the sde file
        with open(output_filename_sde, 'w') as output_file:
            output_file.write(result)

        self.log('Created Sentaurus input file \'{0}\'.'.format(output_filename_sde))

        # *************** Generate the sdevice file based on the template *******************
        # Prepare the substitutions
        node_name_msh = '{0}_msh'.format(os.path.join(output_folder, mesh_name))
        if self._illuminated:
            out_mode = 'light'
        else:
            out_mode = 'dark'
        node_name_out = '{0}_{1}_des'.format(os.path.join(output_folder, mesh_name), out_mode)

        # Change the parfile
        parfile = '{0}_t{1:d}.par'.format(os.path.join(output_folder, 'sdevice'), int(time))
        S0 = self.estimate_srv(cNa=concentration, use_srv=use_srv)

        # Limit to 5 significant digits
        s0_val = format(S0, '1.4e')

        # Prepare the conductivity models for the different shunt regions
        resistivity_template = "Region = \"${shunt_region}\" {\n" \
                               "    Resistivity {\n" \
                               "        * Resist(T) = Resist0 * ( 1 + TempCoef * ( T - 273 ) )\n" \
                               "        Resist0 = ${R0}\n" \
                               "        TempCoef = 1.0E-5\n" \
                               "    }\n" \
                               "}\n\n"
        src = Template(resistivity_template)
        shunt_conductivity = ""
        if dy > 0:
            for i, s in enumerate(sigma):
                if s > 0:
                    r0 = 1.0 / s
                else:
                    r0 = 1E50
                    self.log('Found conductivity = 0 at time {0:1.3f} hr'.format(time / 3600), level='WARNING')
                substitutions = {'shunt_region': 'shunt.Region.{0:d}'.format(i), 'R0': '{0:1.4E}'.format(r0)}
                shunt_conductivity += src.safe_substitute(substitutions)

        substitutions = {
            'S0_val': s0_val,
            'shunt_conductivity': shunt_conductivity,
            'metal_work_function': metal_work_function
        }
        # Load the template file
        template_file_par_file = open(os.path.join(self._cwd, self._parfile_template), 'r')
        src = Template(template_file_par_file.read())
        # perform the substitutions
        result = src.substitute(substitutions)
        template_file_par_file.close()
        # Save the sde file
        with open(parfile, 'w') as output_file:
            output_file.write(result)

        self.log('Created parfile \'{0}\' with SRV value S0_val: {1} and shunt resistivities'.format(parfile, s0_val))

        # Load the template file
        substitutions_sdevice = {
            'nodenum': node_name_msh,
            'parfile': parfile,
            'node_out': node_name_out,
            'illumination_spectrum': os.path.join(self._cwd, self._spectrum),  # Not used, consider removing
            'area_factor': '{0:.3e}'.format(area_factor),
            'folder_path': self._folder_path,
        }
        if three_dimensional:
            template_file_sdevice = open(os.path.join(self._cwd, self._sdevice_template_3d), 'r')
        else:
            substitutions_sdevice['cylindrical'] = symmetry
            template_file_sdevice = open(os.path.join(self._cwd, self._sdevice_template_2d), 'r')
        src = Template(template_file_sdevice.read())
        # perform the substitutions
        result = src.substitute(substitutions_sdevice)
        template_file_sdevice.close()
        output_filename_sdevice = 'sdevice_des_t{0}.cmd'.format(int(time))
        output_filename_sdevice = os.path.join(self._folder_path, output_filename_sdevice)
        # Save the sde file
        with open(output_filename_sdevice, 'w') as output_file:
            output_file.write(result)

        self.log('Created Sentaurus input file \'{0}\'.'.format(output_filename_sdevice))

        return [output_filename_sde, output_filename_sdevice]

    def update_files_traps(self, concentration: np.ndarray, x: np.ndarray,
                           time: Union[int, float], use_srv: bool = False,
                           two_dimensional: bool = True,
                           cylindrical: float = False,
                           activated_na: float = 1.0,
                           fixed_conductivity: float = 0,
                           shunt_depth: float = 0,
                           metal_work_function: float = 4.05) -> List[str]:
        """
        This method modifies the Sentaurus input files to include updated shunt depth and modified external files
        (including the .par file). Creates a shunt with spatially varying conductivity_model based on doping profiles
        calculated at each time.

        Parameters
        ----------
        concentration: np.ndarray
            Sodium concentration from the h5 file
        x: np.ndarray
            depth from the h5 file (list)
        time: Union[int, float]
            time from the h5 file (int or float)
        use_srv: bool
            If true updates the parameter file with the estimated SRV
        two_dimensional: float
            True if a 2D simulation otherwise 3D
        cylindrical: bool
            If True, activate cylindrical geometry
        activated_na: float
            The fraction of sodium that is ionized.
        fixed_conductivity: float
            If the value is 0, ignore it. If > 0 use the value as a constant conductivity in the shunt (Ohm cm)^{-1}
        shunt_depth: float
            If the value is 0, ignore. If >0 use it as shunt depth (in um)
        metal_work_function: float
            The work function of the metal shunt (default 4.05 eV)

        Returns
        -------
        List[str]
            [sde_file_path, sdevice_file_path]
        """
        # The directory to save the data to
        output_folder: str = self._folder_path
        # name for the mesh file
        mesh_name = "n_t{0:d}".format(int(time))
        # Calculate segregation coefficient at each depth and the depth of the shunt
        cseg = np.abs(concentration) * self._segregation_coefficient

        # Get the conductivity
        if fixed_conductivity == 0 and shunt_depth == 0:
            conductivity = self.conductivity_model(cseg, activated_na=activated_na)
            idx_threshold = conductivity >= 1E-3
            conductivity = conductivity[idx_threshold]
            x = x[idx_threshold]

            if time == 0:
                shunt_depth = 0
                dy = 0
                y = np.zeros(1)
                sigma = y
            elif len(x) == 1:
                shunt_depth = x[0]
                dy = shunt_depth
                y = np.zeros(1)
                sigma = y
            else:
                shunt_depth = np.amax(x)  # Depth of the profile in um. The depth does not start at 0 so subtract x[0]

            if len(x) > 100:
                # Interpolate the profile to get only 100 points
                y = np.linspace(0, shunt_depth, 101)
                dy = y[1] - y[0]
                f = interpolate.interp1d(x, conductivity, kind='slinear', fill_value="extrapolate")
                y = y[0:-1]
                sigma = f(y)
            elif len(x) >= 1:
                y = x[0:-1]
                dy = y[1] - y[0]
                sigma = conductivity[0:-1]
            else:
                y = np.zeros(1)
                dy = 0
                sigma = y
        else:
            if time == 0:
                shunt_depth = 0
                dy = 0
                y = np.zeros(1)
                sigma = y
            else:
                shunt_depth = abs(shunt_depth)
                fixed_conductivity = abs(fixed_conductivity)
                y = np.array([0])
                dy = shunt_depth
                sigma = fixed_conductivity * np.ones_like(y)

        box_indices = ' '.join(['%d' % i for i in range(len(y))])

        # *************** Generate the sde file based on the template *******************
        # Prepare the substitutions
        # conductivity_filename = "na_profile_t{0:.0f}.plx".format(time)
        # conductivity_fullname = os.path.join(self._folder_path, conductivity_filename)
        mesh_node_file = os.path.join(self._folder_path, mesh_name)
        substitutions_sde = {
            'dshunt_name': self._dshuntname,
            'shunt_depth': shunt_depth,
            'optical_generation_file': self._optical_generation_file,
            'nodenum': mesh_node_file,
            'contact_length': self._contact_length,
            'device_length': self._cell_width,
            'shunt_positions': ' '.join(map(str, y)),
            'shunt_sigmas': ' '.join(map(str, sigma)),
            'shunt_index': box_indices,
            'shunt_dy': dy
        }
        # Load the template file
        symmetry: str = '* No symmetry'
        if two_dimensional:
            template_file_sde = open(os.path.join(self._cwd, self._sde_template_2d), 'r')
            area_factor = 1E11 / self._cell_width
        else:
            template_file_sde = open(os.path.join(self._cwd, self._sde_template_3d), 'r')
            area_factor = 1E11 / (25 * self._cell_width)

        if cylindrical:
            template_file_sde = open(os.path.join(self._cwd, self._sde_template_3_cylindrical), 'r')
            area_factor = 1E11 / (np.pi * (self._cell_width / 2) ** 2)
            symmetry: str = 'Cylindrical'

        src = Template(template_file_sde.read())
        # perform the substitutions
        result = src.substitute(substitutions_sde)
        template_file_sde.close()
        output_filename_sde = 'sde_dvs_t{0:d}.cmd'.format(int(time))
        output_filename_sde = os.path.join(self._folder_path, output_filename_sde)
        # Save the sde file
        with open(output_filename_sde, 'w') as output_file:
            output_file.write(result)

        self.log('Created Sentaurus input file \'{0}\'.'.format(output_filename_sde))

        # *************** Generate the sdevice file based on the template *******************
        # Prepare the substitutions
        node_name_msh = '{0}_msh'.format(os.path.join(output_folder, mesh_name))
        if self._illuminated:
            out_mode = 'light'
        else:
            out_mode = 'dark'
        node_name_out = '{0}_{1}_des'.format(os.path.join(output_folder, mesh_name), out_mode)
        parfile = '{0}_t{1}.par'.format(os.path.join(output_folder, 'sdevice'), int(time))

        substitutions_sdevice = {
            'nodenum': node_name_msh,
            'parfile': parfile,
            'node_out': node_name_out,
            'illumination_spectrum': os.path.join(self._cwd, self._spectrum),
            'area_factor': '{0:.3e}'.format(area_factor),
            'folder_path': self._folder_path,
            'dshunt_name': self._dshuntname,
            'cylindrical': symmetry
        }

        # Change the parfile

        parfile = '{0}_t{1:d}.par'.format(os.path.join(output_folder, 'sdevice'), int(time))
        S0 = self.estimate_srv(cNa=concentration, use_srv=use_srv)
        # The directory to save the data to
        # output_folder = self._folder_path
        # Limit to 5 significant digits
        s0_val = format(S0, '1.4e')

        # Prepare the conductivity models for the different shunt regions
        resistivity_template = "Region = \"${shunt_region}\" {\n" \
                               "    Resistivity {\n" \
                               "        * Resist(T) = Resist0 * ( 1 + TempCoef * ( T - 273 ) )\n" \
                               "        Resist0 = ${R0}\n" \
                               "        TempCoef = 1.0E-5\n" \
                               "    }\n" \
                               "}\n\n"
        src = Template(resistivity_template)
        shunt_conductivity = ""
        if dy > 0:
            for i, s in enumerate(sigma):
                if s > 0:
                    r0 = 1.0 / s
                else:
                    r0 = 1E50
                    self.log('Found s = 0 at time {0:1.3f} hr'.format(time / 3600), level='WARNING')
                substitutions = {'shunt_region': 'shunt.Region.{0:d}'.format(i), 'R0': '{0:1.4E}'.format(r0)}
                shunt_conductivity += src.safe_substitute(substitutions)

        substitutions = {
            'S0_val': s0_val,
            'shunt_conductivity': shunt_conductivity,
            'metal_work_function': metal_work_function
        }
        # Load the template file
        template_file_par_file = open(os.path.join(self._cwd, self._parfile_template), 'r')
        src = Template(template_file_par_file.read())
        # perform the substitutions
        result = src.substitute(substitutions)
        template_file_par_file.close()
        # Save the sde file
        with open(parfile, 'w') as output_file:
            output_file.write(result)

        self.log('Created parfile \'{0}\' with SRV value S0_val: {1} and shunt resistivities'.format(parfile, s0_val))

        # Load the template file
        template_file_sdevice = open(os.path.join(self._cwd, self._sdevice_template_2d), 'r')
        src = Template(template_file_sdevice.read())
        # perform the substitutions
        result = src.substitute(substitutions_sdevice)
        template_file_sdevice.close()
        output_filename_sdevice = 'sdevice_des_t{0}.cmd'.format(int(time))
        output_filename_sdevice = os.path.join(self._folder_path, output_filename_sdevice)
        # Save the sde file
        with open(output_filename_sdevice, 'w') as output_file:
            output_file.write(result)

        self.log('Created Sentaurus input file \'{0}\'.'.format(output_filename_sdevice))

        return [output_filename_sde, output_filename_sdevice]

    def estimate_srv(self, cNa: Union[np.ndarray, List[float]], Ndop: float = 1E19,
                     use_srv: bool = True, use_na: bool = False) -> float:
        """
        Function giving the parameterization of surface recombination velocity as a function of the surface Na
        concentration

        Parameters
        ----------
        cNa: List[float]
            List containing Na concentration as a function of depth (cm-3)
        Ndop: float
            The doping level at the SiNx/emitter interface
        use_srv: bool
            if  False, S0 is set to zero
        use_na: bool
            If True adds the SRV estimated for Na from experimental lifetime data

        Returns
        -------
        str
            The path to the generated sdevice_tX.par
        """

        if not use_srv:
            return 0  # cm/s
        else:
            # Estimate SRV for p++ emitter
            # Altermatt's model (Eq 7)
            # S_p0 = Sp1*(Ndop/Np1)^(gp1) + Sp2*(Ndop / Np2)^(gp2)
            # Sp1 = 500, Sp2 = 60 (cm/s) # Untextured
            # Sp1 = 2800, Sp2 = 300 (cm/s) # Textured
            # gp1 = 0.6, gp2 = 3
            # Np1 = Np2 = 1E19 cm^-3
            gamma1, gamma2 = 0.6, 3
            Sp1, Sp2 = 2800, 300
            Np1, Np2 = 1E19, 1E1
            Sp0 = Sp1 * (Ndop / Np1) ** gamma1 + Sp2 * (Ndop / Np1) ** gamma2

            # Limit S0 to the thermal velocity of electrons

            me = 9.1e-31  # electron mass, kg
            kB = 1.38e-23  # J.K-1
            TK = self._temperature + 273.15

            # (cm/s) thermal energy for non-relastivistic electrons E=df*kB*T where df number of degrees of freedom,
            # and E=m*v^2
            vth = 100 * np.sqrt(3 * kB * TK / me)

            if use_na:
                # Calculate surface recombination velocity
                # Fitting values according to the phosphorus parameterization by
                # Altermatt et al, Journal of App Phys 92, 3187, 2002.
                S1 = 500  # cm/s
                S2 = 60  # cm/s
                N1 = 1e10  # cm-3 (modified from Altermatt et al)
                N2 = 1e10  # cm-3 (modified from Altermatt et al)

                # Parameterization of the surface recombination velocity
                # TODO: verify how to implement this. The values N1, N2 seem to be coming from life time data from ASU
                S0_Na = S1 * (cNa[0] / N1) ** gamma1 + S2 * (
                        cNa[0] / N2) ** gamma2  # Altermatt et al, Journal of App Phys 92, 3187, 2002
            else:
                S0_Na = 0

            S0 = Sp0 + S0_Na

            if S0 > vth:  # cm/s
                S0 = vth  # cm/s

            return S0

    def srv_param(self, cNa: Union[np.ndarray, List[float]], time: Union[int, float],
                  use_srv: bool) -> str:
        """
        Function giving the parameterization of surface recombination velocity as a function of the surface Na
        concentration

        Parameters
        ----------
        cNa: List[float]
            List containing Na concentration as a function of depth (cm-3)
        time: Union[int, float]
            Temperature (C) (int or float)
        use_srv: bool
            if  False, S0 is set to zero

        Returns
        -------
        str
            The path to the generated sdevice_tX.par
        """

        S0 = self.estimate_srv(cNa=cNa, use_srv=use_srv)

        # The directory to save the data to
        output_folder = self._folder_path
        # Limit to 5 significant digits
        s0_val = format(S0, '1.4e')
        parfile = '{0}_t{1:d}.par'.format(os.path.join(output_folder, 'sdevice'), int(time))
        substitutions = {
            'S0_val': s0_val,
        }
        # Load the template file
        template_file_par_file = open(os.path.join(self._cwd, self._parfile_template), 'r')
        src = Template(template_file_par_file.read())
        # perform the substitutions
        result = src.substitute(substitutions)
        template_file_par_file.close()
        # Save the sde file
        with open(parfile, 'w') as output_file:
            output_file.write(result)

        self.log('Created parfile \'{0}\' with SRV value S0_val: {1}.'.format(parfile, s0_val))
        return parfile

    def read_jv(self, h5_filename: str) -> np.ndarray:
        """
        Reads the IV curve from the plt file generated by sentaurus

        Parameters
        ----------
        h5_filename: str
            The basename .plt JV file

        Returns
        -------
        np.ndarray
            An array with the current density as a function of the voltage
        """
        full_file = os.path.join(self._folder_path, h5_filename)
        # main_group = 'collection/geometry_0/state_0'
        voltage_dataset_name = 'em_contact OuterVoltage'
        current_dataset_name = 'base_contact TotalCurrent'
        voltage_dataset = self.tdr_get_plt_dataset(h5file=full_file, dataset_name=voltage_dataset_name)
        current_dataset = self.tdr_get_plt_dataset(h5file=full_file, dataset_name=current_dataset_name)

        voltage = np.array(voltage_dataset)
        current = np.array(current_dataset)
        jv = np.array(
            [(v, j) for v, j in zip(voltage, current)]
        )
        return jv

    def plt2tdr(self):
        """
        Converts the output the DF-ISE plt files to tdr (HDFS) files
        """
        data_folder = self._folder_path
        output_folder = os.path.join(data_folder, 'jv_plots')
        file_index_csv = 'file_index.csv'
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        file_index = np.empty(len(self._simulation_time), dtype=np.dtype([
            ('time (s)', 'i'), ('filename', 'U100'), ('index', 'i')
        ]))

        for i, t in enumerate(self._simulation_time):
            file_name = 'n_t{0:d}_light_des.plt'.format(int(t))
            output_file = 'n_t{0:d}_light_des.tdr'.format(int(t))
            subprocess.run(['tdx', '-d', os.path.abspath(os.path.join(data_folder, file_name)),
                            os.path.join(output_folder, output_file)])
            file_index[i] = (t, output_file, self._simulation_indices[i])

        df = pd.DataFrame(data=file_index)
        df.to_csv(path_or_buf=os.path.join(output_folder, file_index_csv), index=False)

    def export_xcut_band_diagram(self, tdr_file, x_cut: float):
        """
        Creates a cut along the y axis of the tdr mesh and exports the band diagram

        Parameters
        ----------
        tdr_file: str
            The path to the 2D file with the datasets
        x_cut: float
            The x-position of the cut along the y axis in um

        Returns
        -------

        """
        datasets = ['Y', 'ConductionBandEnergy', 'ValenceBandEnergy', 'eQuasiFermiEnergy', ' hQuasiFermiEnergy']
        self.export_x_cut(tdr_file=tdr_file, x_cut=x_cut, variables=datasets)

    def export_x_cut(self, tdr_file: str, x_cut: float, variables: List[str]):
        """
        Creates a cut along the y axis of the tdr mesh and exports the variables in the list

        Parameters
        ----------
        tdr_file: str
            The path to the 2D tdr file with the datasets
        x_cut: float
            The x-position of the cut along the y axis in um
        variables: List[str]
            The list of variables to export

        Returns
        -------

        """
        # The directory to save the data to
        output_folder = self._folder_path
        base_filename = os.path.splitext(os.path.basename(tdr_file))[0]

        batch_file = '{0}_{1}.tcl'.format(os.path.join(output_folder, 'batch_bd'), base_filename)

        substitutions = {
            'tdr_file': tdr_file,
            'x_cut_line': x_cut,
            'variables': ' '.join(map(str, variables))
        }
        # Load the template file
        template_file_par_file = open(os.path.join(self._cwd, self._batch_band_diagram_template), 'r')
        src = Template(template_file_par_file.read())
        # perform the substitutions
        result = src.substitute(substitutions)
        template_file_par_file.close()
        # Save the sde file
        with open(batch_file, 'w') as output_file:
            output_file.write(result)

        self.log('Created svisual batch process \'{0}\'.'.format(batch_file))

        os.system('svisual -batch {0}'.format(batch_file))

    @property
    def na_profile_time_array(self):
        return self._na_profile_time_array

    @staticmethod
    def tdr_list_plt_datasets(h5file: str) -> dict:
        main_group = 'collection/geometry_0/state_0'
        with h5py.File(h5file, 'r') as hf:
            hf_datasets = list(hf[main_group])
            ds = {hf['collection/geometry_0/state_0'][d].attrs['name'].astype(str): d for d in hf_datasets}
        return ds

    def tdr_get_plt_dataset(self, h5file: str, dataset_name: str):
        available_datasets = self.tdr_list_plt_datasets(h5file=h5file)
        main_group = 'collection/geometry_0/state_0'
        if dataset_name in available_datasets:
            with h5py.File(h5file, 'r') as hf:
                ds = hf[main_group][dataset_name]
        else:
            ds = None
        return ds

    def run_sde(self, input_file: str) -> bool:
        try:
            os.system('sde -e -l {0}'.format(input_file))  # subprocess.run(['sde', '-e', '-l', input_file], shell=True)
        except FileNotFoundError as e:
            self.log(e.strerror)
            return False
        else:
            return True

    def run_sdevice(self, input_file: str) -> bool:
        try:
            os.system('sdevice {0}'.format(input_file))
            # r = subprocess.run(['sdevice', input_file], shell=True)
            # if r.returncode != 0:
            #     self.log('Error running sdevice.')
            #     return False
        except FileNotFoundError as e:
            self.log(e.strerror)
            return False
        else:
            return True

    @staticmethod
    def run_svisual(mesh_file: str):
        subprocess.Popen(['svisual', mesh_file], shell=False)

    def check_errfile(self) -> bool:
        """
        Checks errfile to see whether sde execution caused an error

        Returns
        -------
        bool
            True if there is an error, False otherwise
        """
        error_file = os.path.join(os.path.dirname(self._cwd), 'errfile.txt')
        with open(error_file, 'r') as f:
            error_flag = f.read()
        return not bool(error_flag)

    def __create_logger(self, logger_id: str = 'pid_logger') -> logging.Logger:
        logger_filename = os.path.join(self._folder_path, 'pidlogger.log')
        pid_logger = logging.getLogger(logger_id)
        pid_logger.setLevel(logging.DEBUG)
        # create file handler which logs even debug messages
        fh = logging.FileHandler(logger_filename)
        fh.setLevel(logging.DEBUG)
        # create console handler and set level to debug
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        # create formatter and add it to the handlers
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        fh.setFormatter(formatter)

        # Make sure no other handlers exist
        pid_logger.handlers = []
        # add the handlers to the logger
        pid_logger.addHandler(fh)  # < handlers[0] (In case the logfile needs to be retrived)
        pid_logger.addHandler(ch)  # < handlers[1]
        return pid_logger

    def prepare_folder(self, folder: str, overwrite: bool = True):
        """
        Prepare the folder structure. 
        
        It will check if the folder exists before attempting to create it.
        If the the overwrite option is True it will not attempt to create a new folder.
        Otherwise. it will append a consecutive number to the folder name.
        
        Parameters
        ----------
        folder: str
            The name of the folder to be created. 
        overwrite: bool
            If set to False it will append a consecutive number to the folder name.

        Returns
        -------

        """
        if not overwrite:
            folder_path = utils.make_new_folder(folder)
        else:
            folder_path = folder
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
        opt_gen_file = os.path.basename(self._optical_generation_file)

        shutil.copy(self._optical_generation_file, os.path.abspath(os.path.join(folder_path, opt_gen_file)))
        return folder_path

    def log(self, msg: str, level="INFO"):
        """

        Parameters
        ----------
        msg: str
            The message to print
        level: level
            The level of the message
        """
        level_no = self._loggingLevels[level]
        if self._logger is None:
            print(msg)
        else:
            self._logger.log(level_no, msg)
