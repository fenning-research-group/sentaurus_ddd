#!/usr/bin/env python3
# from typing import Union, List
import pid.pidmodel as pdm
import pid.utils as utils
import numpy as np
import argparse
import configparser
# from scipy import special
import h5py
import os

pi_sqrt = np.sqrt(np.pi)
emitter_depth = 0.3  # um
# c_min_at_junction = 1E16  # cm^-3
sigma_min_junction = 1

# def f(t: Union[np.ndarray, List[float], float], x: float, target_c: float,
#       c_s: float, diffusivity: float) -> Union[np.ndarray, float]:
#     return c_s * special.erfc(x / (2.0 * np.sqrt(diffusivity * t))) - target_c
#
#
# def fp(t: Union[np.ndarray, List[float], float], x: Union[np.ndarray, List[float], float],
#        c_s: float, diffusivity: float) -> Union[np.ndarray, float]:
#     a = x / np.sqrt(2.0 * diffusivity)
#     return c_s * a * np.exp(-a ** 2 / t) * t ** (-3 / 2) / pi_sqrt
#
#
# def solve_diffusion_at_x(c_s: float, diffusivity: float, position: float, concentration: float = 1E16, **kwargs):
#     iterations: int = kwargs.get('iterations', 500)
#     w: float = kwargs.get('damping', 1.0)
#     t = 0.1
#     for i in range(iterations):
#         t = t - w * f(t, position, concentration, c_s, diffusivity) / fp(t, position, c_s, diffusivity)
#     return t


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='The input file.',
                        required=True)
    parser.add_argument('-u', '--unattended', help='Unattended execution.',
                        action='store_true',
                        required=False)
    parser.add_argument('-s', '--source_limited', help='Source limited simulation.',
                        action='store_true',
                        required=False)
    args = parser.parse_args()
    config_file = args.input
    attended = not args.unattended
    source_limited = args.source_limited


    # Instantiate the config parser
    config = configparser.ConfigParser()

    print('Reading configuration file...')
    # Load the configuration file
    config.read(config_file)

    h5file = config.get(section='transport', option='hf_file')
    sc = config.getfloat(section='transport', option='segregation_coefficient')

    cell_width = config.getfloat(section='device', option='cell_width')
    cell_length = config.getfloat(section='device', option='cell_length')
    use_srv = config.getboolean(section='device', option='use_srv')
    time_sampling = config.get(section='device', option='time_sampling')
    time_steps = config.getint(section='device', option='time_steps')
    t_max = config.getfloat(section='device', option='simulation_time')
    shunt_depth = config.getfloat(section='device', option='shunt_depth', fallback=0)
    simulate_zero = config.getboolean(section='device', option='simulate_t0', fallback=True)
    conductivity_cutoff = config.getfloat(section='device', option='conductivity_cutoff', fallback=1E-10)
    print('Simulating t = 0: {0}'.format(simulate_zero))

    rho = 1. / (cell_width * cell_length)

    with h5py.File(h5file, 'r') as hf:
        dst = hf['/time']
        if source_limited: # spyder
            n_na = dst.attrs['na_0']
        else:
            cs = dst.attrs['Csource']
        d1 = hf['/L1'].attrs['D']
        d2 = hf['/L2'].attrs['D']
        mtc = dst.attrs['h']
        t_sim = np.array(hf.get('time'))
        depth = np.array(hf.get('L2/x'))

    simulation_thickness = np.amax(depth - depth[0])
    x_max = np.amax(emitter_depth)
    dt = t_sim[2] = t_sim[1]
    # find the index in the depth array corresponding to (approximately) the emitter thickness
    idx_emitter = (np.abs(depth - emitter_depth)).argmin()

    if shunt_depth == 0:
        shunt_depth = simulation_thickness * np.sqrt(2)

    shunt_length = shunt_depth * np.sqrt(2)
    print('Shunt length = {0:.3f} um.'.format(shunt_length))
    print('Shunt depth = {0:.3f} um.'.format(shunt_depth))

    results_path = '/home/erickmtz/sentaurus_pid/results/3D/recovery/'
    if source_limited:
        results_path += 'N_Na={0:.1E}moles_D1=4E-16_h={1:.0E}_D2={2:.0E}_rho={3:.0E}_SD={4:.1f}_s={5:.0E}'.format(
            n_na, mtc, d2, rho, shunt_depth, sc
        )
    else:
        results_path += 'Cs={0:.0E}D1=4E-16_h={1:.0E}_D2={2:.0E}_rho={3:.0E}_SD={4:.1f}_s={5:.0E}'.format(
            cs, mtc, d2, rho, shunt_depth, sc
        )
    print('Results path:')
    print(results_path)
    pid_model = pdm.PIDModel(
        folder_path=results_path, na_h5_file=h5file, segregation_coefficient=sc, use_srv=use_srv,
        cell_width=cell_width, cell_length=cell_length
    )

    # Estimate the the time at which a concentration of 10E16 will reach the emitter using a simple erf model
    t_shunting = 0
    dc = np.inf
    for i, c in enumerate(t_sim):
        ct_ds = 'ct_{0:d}'.format(i)
        with h5py.File(h5file, 'r') as hf:
            group_si_c = hf['/L2/concentration']
            ct = np.abs(np.array(group_si_c[ct_ds]))
            # c_min = np.amin(np.abs(ct))
            # ct[ct < 0] = c_min
            # if np.any(ct < 0):
            #     print("There are negative values of the concentration!")
            #     print(ct)
            sigma = utils.conductivity_model(concentration=ct[idx_emitter], segregation_coefficient=sc)

            d = abs(sigma - sigma_min_junction)
            if d < dc:
                dc = d
                t_shunting = t_sim[i]

    t_min = 0 if simulate_zero else dt
    print('Profile time series dt: {0:.3f} s'.format(dt))
    print('Starting sumulation at time: {0:.1f} s'.format(t_min))
    print('Maximum simulation time: {0:.1f} h'.format(t_max / 3600))
    print('Time to reach {0:.1E} S/cm at the junction: {1:.1f} s'.format(sigma_min_junction, t_shunting))
    # If the time required to form the shunt is less than the simulation time
    # Sample the time points until shunt formation faster and slower beyond that point

    requested_time_points = pid_model.request_time_points(
        t_min=0, t_max=t_max, t_steps=time_steps, spacing=time_sampling
    )

    req_sim_times = np.array([t_sim[rt] / 3600 for rt in requested_time_points])

    print('Simulation times: ')
    print(req_sim_times)
    print('Simulation time indices: ')
    print(requested_time_points)

    if attended:
        input("Press Enter to continue...")

    pid_model.run_pid_metal(
        requested_indices=requested_time_points, overwrite_folder=True, use_srv=use_srv,
        conductivity_cutoff=conductivity_cutoff,
        shunt_geometry='triangle', shunt_depth=shunt_depth, n_regions=int(shunt_length/0.01)
    )

    # Get the basename of the config file
    config_file_name = os.path.basename(config_file)
    # Construct the path to the working directory
    conifg_file_destination = os.path.join(results_path), config_file_name
    # Move the input file to the working folder
    os.system('mv {0} {1}'.format(config_file, conifg_file_destination))

    # Copy the results to the host computer (remember to delete them later)
    print('Copying results to the host file system...')
    guest_folder = 'results/3D/*'
    host_folder = r'/home/erickmtz/shared/centos7/sentaurus_pid/results/3D/'
    guest_key = '/home/erickmtz/.ssh/erick_guest_rsa'
    os.system('scp -r -i {guest_key} {guest_folder} erickmtz@mae202-34.ucsd.edu:{host_folder}'.format(
        guest_key=guest_key,
        guest_folder=guest_folder,
        host_folder=host_folder
    ))

