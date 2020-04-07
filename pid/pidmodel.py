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

na_profile_type = np.dtype([('x (um)', 'd'), ('c (1/cm3)', 'd')])
jv_dtype = np.dtype([('voltage (V)', 'd'), ('current (mA/cm2)', 'd')])
na_file = "two_layers_D1=4E-16cm2ps_D2=1E-14cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-10_m1.0e+00_pnp.h5"

class PIDModel:
    """
    This class provides methods to run the electrical module of the DDD model
    """
    _NaProfileFileName: str = None
    _cSi: float = 5e22
    _dshuntname: str = "dshunt1"
    _logger: logging.Logger = None
    _loggingLevels: dict = {'NOTSET': logging.NOTSET,
                            'DEBUG': logging.DEBUG,
                            'INFO': logging.INFO,
                            'WARNING': logging.WARNING,
                            'ERROR': logging.ERROR,
                            'CRITICAL': logging.CRITICAL}
    _sde_template: str = 'sde_dvs.scm'
    _sdevice_template_dark: str = 'sdevice_dark_des.cmd'
    _sdevice_template_light: str = 'sdevice_light_des.cmd'
    _simulation_time: np.ndarray = np.array([])
    _simulation_indices: np.ndarray = np.array([])
    _spectrum: str = 'am15-g_reduced.txt'
    _parfile_template: str = 'sdevice.txt'

    def __init__(self, folder_path: str = "/home/linux/ieng6/na299x/na299x/DB/Erick/PID/test",
                 na_h5_file: str = na_file,
                 segregation_coefficient: float = 10,
                 use_srv: bool = False,
                 clathrate_file: str = "clathrates_cond.csv",
                 illuminated: bool = True,
                 contact_length: float = 5,
                 device_length: float = 5):
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
        illuminated: bool
            True if the device is under illumination
            TODO: Consider that the device is always under illumination
        contact_length: float
            The length of the metal contact
        device_length: float
            The total length of the cell
        """

        self._clathrateFile: str = clathrate_file
        self._cwd: str = os.path.dirname(__file__)
        self._folder_name: str = folder_path
        self._folder_path: str = folder_path
        self._NaProfileFileName: str = na_h5_file
        self._segregation_coefficient: float = segregation_coefficient
        self._useSRV: bool = use_srv
        self._sdevice_template = self._sdevice_template_light if illuminated else self._sdevice_template_dark
        self._illuminated = illuminated
        self._contact_length = contact_length
        self._device_length = device_length
        with h5py.File(na_h5_file, 'r') as hf:
            group_time = hf['time']
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
        self.log('Sentaurus editor template file: {0}'.format(self._sde_template))
        self.log('Sentaurus device template file: {0}'.format(self._sdevice_template))

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
            conductivity = self.conductivity_model(ct[:])
            # Save conductivity profile in a .plx file
            conductivity_filename = "conductivity_t{0:d}.plx".format(time[i])
            conductivity_fullname = os.path.join(self._folder_path, conductivity_filename)
            with open(conductivity_fullname) as fp:
                fp.write("\"MetalConductivity\"\n")
                for xi, c in zip(x, conductivity):
                    line = "{0:1.4E}\t{1:1.4E}\n".format(xi, c)
                    fp.write(line)
            self.log('Created concentration file \'{0}\'.'.format(conductivity_fullname))
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

    def conductivity_model(self, concentration: np.ndarray):
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

        Returns
        -------
        np.ndarray
            The conductivity_model profile
        """

        # Na concentration in the shunt
        cshunt = concentration * self._segregation_coefficient

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
        This method modifyies the Sentaurus input files to include updated shunt depth and modified external files
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
        mesh_name = "n_t{0:d}".format(time)
        # Calculate segregation coefficient at each depth and the depth of the shunt
        cseg = concentration * self._segregation_coefficient
        L = len(cseg)
        # pdb.set_trace()
        k = 0
        # ATTENTION: need to modify this as the Na concentration is not a meaningful value here. The shunt depth is
        # highly sensitive to the conductivity_model model and to segregation_coefficient.
        while cseg[k] > 0.1 and k < L:
            # Consider only the part of the [Na] profile in the stacking fault that is at least 1 cm-3
            # (to remove useless extra points). TO CHANGE.
            k += 1

        # **************************************************************************************************************
        # ATTENTION! This should be much more reliable than using cseg[k]>0.1, to be tested.
        # Or maybe we shouldn't limit the shunt depth at all.
        # while(PIDmodel.conmodel(cseg[k],60,mseg) > 1e-3 and k<L-1): # Limit the shunt depth to conductivities higher
        # than 1e-3 S/cm (assuming 1e-3 is lower than the limit conductivity_model for shunting)
        # k=k+1
        # **************************************************************************************************************
        shunt_depth = x[k] - x[0]  # Depth of the profile in um. The depth does not start at 0 so subtract x[0]

        # Define an arbitrary shunt depth in case of very low concentrations
        # that would lead to a shunt depth of 0 um and cause Sentaurus to crash
        if abs(shunt_depth) < 1E-4:
            shunt_depth = 0.1

        # NOTE: CHECK IF IT MATTERS IF THE EXTERNAL PROFILE IS DEEPER THAN THE DEPTH DEFINED HERE. CURRENTLY IT IS.
        # *************** Generate the sde file based on the template *******************
        # Prepare the substitutions
        conductivity_filename = "conductivity_t{0:d}.plx".format(time)
        conductivity_fullname = os.path.join(self._folder_path, conductivity_filename)
        mesh_node_file = os.path.join(self._folder_path, mesh_name)
        substitutions_sde = {
            'dshunt_name': self._dshuntname,
            'shunt_depth': shunt_depth,
            'shunt_conductivity_file': conductivity_fullname,
            'nodnum1': mesh_node_file,
            'contact_length': self._contact_length,
            'device_length': self._device_length
        }
        # Load the template file
        template_file_sde = open(os.path.join(self._cwd, self._sde_template), 'r')
        src = Template(template_file_sde.read())
        # perform the substitutions
        result = src.substitute(substitutions_sde)
        template_file_sde.close()
        output_filename_sde = 'sde_dvs_t{0}.cmd'.format(int(time))
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
            'nodenum1': node_name_msh,
            'parfile': parfile,
            'nodeout': node_name_out,
            'illumination_spectrum': os.path.join(self._cwd, self._spectrum),
            'area_factor': '{0:.3e}'.format(1E11/self._device_length)
        }
        # Load the template file
        template_file_sdevice = open(os.path.join(self._cwd, self._sdevice_template), 'r')
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

        if not use_srv:
            S0 = 0  # cm/s
        else:
            # Calculate surface recombination velocity
            # Fitting values according to the phosphorus parameterization by
            # Altermatt et al, Journal of App Phys 92, 3187, 2002.
            S1 = 500  # cm/s
            S2 = 60  # cm/s
            gamma1 = 0.6
            gamma2 = 3
            N1 = 1e10  # cm-3 (modified from Altermatt et al)
            N2 = 1e10  # cm-3 (modified from Altermatt et al)

            # Parameterization of the surface recombination velocity
            S0 = S1 * (cNa[0] / N1) ** gamma1 + S2 * (
                        cNa[0] / N2) ** gamma2  # Altermatt et al, Journal of App Phys 92, 3187, 2002

            # Limit S0 to the thermal velocity of electrons

            me = 9.1e-31  # electron mass, kg
            kB = 1.38e-23  # J.K-1
            TK = self._temperature + 273.15

            # (cm/s) thermal energy for non-relastivistic electrons E=df*kB*T where df number of degrees of freedom,
            # and E=m*v^2
            vth = 100 * np.sqrt(3 * kB * TK / me)

            if S0 > vth:  # cm/s
                S0 = vth  # cm/s

        # The directory to save the data to
        output_folder = self._folder_path
        # Limit to 5 significant digits
        s0_val = format(S0, '1.4e')
        parfile = '{0}_t{1:d}.par'.format(os.path.join(output_folder, 'sdevice'), time)
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
            ('time (s)', 'i'), ('filename', 'S[100]'), ('index', 'i')
        ]))

        for i, t in enumerate(self._simulation_time):
            file_name = 'n_t{0:d}_light_des.plt'.format(t)
            output_file = 'n_t{0:d}_light_des.tdr'.format(t)
            subprocess.run(['tdx', '-d', os.path.join(data_folder, file_name),
                                os.path.join(output_folder, output_file)])
            file_index[i] = (t, output_file, self._simulation_indices[i])

        df = pd.DataFrame(data=file_index)
        df.to_csv(path_or_buf=os.path.join(output_folder, file_index_csv), index=False)



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
            r = subprocess.run(['sde', '-e', '-l', input_file], shell=True)
            if r.returncode != 0:
                self.log('Error running sde.')
                return False
        except FileNotFoundError as e:
            self.log(e.strerror)
            return False
        else:
            return True

    def run_sdevice(self, input_file: str) -> bool:
        try:
            r = subprocess.run(['sdevice', input_file], shell=True)
            if r.returncode != 0:
                self.log('Error running sdevice.')
                return False
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
        error_file = os.path.join(self._folder_path, 'errfile.txt')
        with open(error_file, 'r') as f:
            error_flag = f.read()

        return bool(error_flag)

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

    @staticmethod
    def prepare_folder(folder: str, overwrite: bool = False):
        if not overwrite:
            folder_path = utils.make_new_folder(folder)
        else:
            folder_path = folder
            os.makedirs(folder_path)
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
