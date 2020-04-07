import os
import h5py
import numpy as np
import pandas as pd
import pid.utils

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import ScalarFormatter
import platform
import h5py
import os
from typing import List
from scipy import interpolate
import matplotlib.gridspec as gridspec


jv_dtype = np.dtype([('voltage (V)', 'd'), ('current (mA/cm2)', 'd')])
nu_t_dtype = np.dtype([
    ('time (s)', 'd'),
    ('jsc (mA/cm2)', 'd'),
    ('voc (V)', 'd'),
    ('pd_mpp (mW/cm2)', 'd'),
    ('v_mpp (V)', 'd'),
    ('j_mpp (mA/cm2)', 'd'),
    ('efficiency', 'd')
])

defaultPlotStyle = {'font.size': 14,
                    'font.family': 'Arial',
                    'font.weight': 'regular',
                    'legend.fontsize': 14,
                    'mathtext.fontset': 'stix',
                    'xtick.direction': 'in',
                    'ytick.direction': 'in',
                    'xtick.major.size': 4.5,
                    'xtick.major.width': 1.75,
                    'ytick.major.size': 4.5,
                    'ytick.major.width': 1.75,
                    'xtick.minor.size': 2.75,
                    'xtick.minor.width': 1.0,
                    'ytick.minor.size': 2.75,
                    'ytick.minor.width': 1.0,
                    'ytick.right': False,
                    'lines.linewidth': 2.5,
                    'lines.markersize': 10,
                    'lines.markeredgewidth': 0.85,
                    'axes.labelpad': 5.0,
                    'axes.labelsize': 16,
                    'axes.labelweight': 'regular',
                    'legend.handletextpad': 0.2,
                    'legend.borderaxespad': 0.2,
                    'axes.linewidth': 1.25,
                    'axes.titlesize': 16,
                    'axes.titleweight': 'bold',
                    'axes.titlepad': 6,
                    'figure.titleweight': 'bold',
                    'figure.dpi': 100}


class Analysis:

    def __init__(self, folder_path: str):
        self._folder_path = folder_path

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

    def batch_efficiency(self, csv_index: str) -> pd.DataFrame:
        """
        Batch method to estimate the efficiency of multiple jv_plots

        Parameters
        ----------
        csv_index: str
            The path to the csv file containing the index of tdr plots linked to each simulation time

        Returns
        -------

        """
        df = pd.read_csv(filepath_or_buffer=csv_index)
        folder_path = os.path.dirname(csv_index)
        results_path = os.path.join(folder_path, 'efficiency_results.csv')
        results = np.empty(len(df), dtype=nu_t_dtype)
        for i, row in df.iterrows():
            t = float(row['time (s)'])
            jv = self.read_jv(h5_filename=row['filename'])
            nu_data = self.estimate_efficiency(voltage=jv['voltage (V)'], current=jv['current (mA/cm2)'])
            results[i] = (
                t, nu_data['jsc'], nu_data['voc'], nu_data['pd_mpp'], nu_data['v_mpp'], nu_data['j_mpp'],
                nu_data['efficiency']
            )

        df_results = pd.DataFrame(data=results)
        df_results.to_csv(path_or_buf=results_path, index=False)
        return df_results

    def plot_jv_t(self, csv_index: str, na_file: str, color_map: str = 'viridis'):
        """
        Plots the JV curves for different times

        Parameters
        ----------
        csv_index: str
            The path to the csv file with the index of tdr plots
        na_file: str
            The path to the h5file with the concentration profiles
        color_map: str
            The name of the color map to use

        Returns
        -------

        """

        efficiency_data = self.batch_efficiency(csv_index=csv_index)
        time_s = efficiency_data['time_s']
        time_h = time_s/3600
        df = pd.read_csv(filepath_or_buffer=csv_index)
        results_path = os.path.dirname(csv_index)

        mpl.rcParams.update(defaultPlotStyle)

        cmap = mpl.cm.get_cmap(color_map)
        normalize = mpl.colors.Normalize(vmin=np.amin(time_h), vmax=np.amax(time_h))
        jv_colors = [cmap(normalize(t)) for t in time_h]
        scalar_maps = mpl.cm.ScalarMappable(cmap=cmap, norm=normalize)

        with h5py.File(na_file, 'r') as hf:
            group_si = hf['L2']
            x = np.array(group_si['x'])
            x = x - x[0]

            # D1 = hf['L1'].attrs['D']
            # D2 = hf['L2'].attrs['D']
            # E = hf['L1'].attrs['electric_field_eff']
            # TC = hf['time'].attrs['temp_c']

        fig = plt.figure()
        fig.set_size_inches(6.5, 3, forward=True)
        fig.subplots_adjust(hspace=0.15, wspace=0.45)
        gs0 = gridspec.GridSpec(ncols=1, nrows=2, figure=fig, width_ratios=[1])
        gs00 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2,
                                                subplot_spec=gs0[0])
        gs10 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1,
                                                subplot_spec=gs0[1])

        ax1 = fig.add_subplot(gs00[0, 0])
        ax2 = fig.add_subplot(gs00[0, 1])
        ax3 = fig.add_subplot(gs10[0, 0])

        with h5py.File(na_file, 'r') as hf:
            group_si = hf['L2']
            profiles = group_si['concentration']
            for i, row in df.iterrows():
                ct_ds = 'ct_{0:d}'.format(row['index'])
                c = np.array(profiles[ct_ds])
                jv = self.read_jv(h5_filename=row['filename'])
                ax1.plot(jv['voltage (V)'], jv['current (mA/cm2)'], '-', color=jv_colors[i])
                ax3.plot(x, c, color=jv_colors[i])

        # interpoate the efficiency
        f = interpolate.interp1d(time_h, efficiency_data['efficiency'], kind='spline')
        t_interp = np.linspace(np.amin(time_h), np.amax(time_h), num=500)
        nu_interp = f(t_interp)
        ax2.scatter(time_h, efficiency_data['efficiency']*1000, marker='o', facecolors='none', label='Simulation')
        ax2.plot(time_h, nu_interp, ':', color='k', dashes=(5, 6), label='Guide to\nthe eye')

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="7.5%", pad=0.03)
        cbar = fig.colorbar(scalar_maps, cax=cax)
        cbar.set_label('Time (h)', rotation=90)

        ax1.set_xlabel('Bias (V)')
        ax1.set_ylabel('J (mA/cm$^2$)')

        ax2.set_xlabel('Time (h)')
        ax2.set_ylabel('Efficiency (%)')

        ax3.set_xlabel('Depth (um)')
        ax3.set_ylabel('[Na] (cm$^{-3}$)')
        ax3.set_yscale('log')

        ax2.legend(loc='lower left', frameon=False)

        ax1.xaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
        ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator(4))
        ax1.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
        ax1.yaxis.set_minor_locator(mticker.AutoMinorLocator(4))

        ax2.xaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
        ax2.xaxis.set_minor_locator(mticker.AutoMinorLocator(4))
        ax2.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
        ax2.yaxis.set_minor_locator(mticker.AutoMinorLocator(4))

        locmaj = mpl.ticker.LogLocator(base=10.0, numticks=5, subs=(1.0, ))
        locmin = mpl.ticker.LogLocator(base=10.0, numticks=50, subs=np.arange(2, 10) * .1)

        ax3.set_ylim(bottom=1E12)
        ax3.yaxis.set_major_locator(locmaj)
        ax3.yaxis.set_minor_locator(locmin)

        fig.savefig(os.path.join(results_path, 'efficiency_t.png'), dpi=600)

        plt.tight_layout()
        plt.show()

    @staticmethod
    def estimate_efficiency(voltage: np.ndarray, current: np.ndarray) -> dict:
        """
        Estimates the efficiency of a cell from its illuminated JV characteristic

        Parameters
        ----------
        voltage: np.ndarray
            The voltage in V
        current: np.ndarray
            The current density in mA/cm2

        Returns
        -------
        dict:
            A dictionary containing the efficiency parameters
        """
        # Estimate the power density
        power_density = current * voltage
        # Find max powerpoint
        pd_mpp = np.amax(power_density)
        idx_mpp = (np.abs(power_density - pd_mpp)).argmin()
        # Find the V_mpp
        v_mpp = voltage[idx_mpp]
        # Find J_mpp
        j_mpp = current[idx_mpp]
        # Find Jsc
        jsc = current[(np.abs(voltage - 0)).argmin()]
        # Find Voc
        voc = voltage[(np.abs(current - 0)).argmin()]
        # Find efficiency
        nu = pd_mpp / 100  # Integrated AM1.5g power density 100 mA/cm2
        # wrap everything in a dictionary
        result = {
            'jsc': jsc,
            'voc': voc,
            'pd_mpp': pd_mpp,
            'v_mpp': v_mpp,
            'j_mpp': j_mpp,
            'efficiency': nu
        }
        return result

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
