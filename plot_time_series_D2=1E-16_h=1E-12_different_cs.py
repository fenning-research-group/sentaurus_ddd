# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 13:15:45 2020

@author: Erick
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import ScalarFormatter
import platform
import os
import matplotlib.gridspec as gridspec
from scipy import interpolate
import h5py
import pid.analysis as pia


base_folder = r'G:\Shared drives\FenningLab2\Projects\PVRD1\Simulations\Sentaurus PID\results\3D'
t_max = 96
simulation_result_dirs = [
        r'Cs=1E+18D1=4E-16_h=1E-12_D2=1E-16_rho=4E-05_SD=1.0_s=5E+01',
        r'Cs=1E+22D1=4E-16_h=1E-12_D2=1E-16_rho=4E-05_SD=1.0_s=5E+01',
    ]

sweep_variable = r'$C_s$ (cm$^{-3}$)'
sweep = [1E18, 1E22]
sweep_log = True

results_folder = r'G:\Shared drives\FenningLab2\Projects\PVRD1\Simulations\Sentaurus PID\results\analysis\0.5MVcm_D2=1E-16_h=1E-12_variable_Cs'

literature_files = [
    {
        'file': r'G:\My Drive\Research\PVRD1\Literature\PID_degradation_time\Masuda2016_Fig5.csv',
        'time_units': 'min',
        'label': 'Masuda 2016',
        'color': 'tab:red',
        'marker': 'o',
        'type': 'power',
        'normalized': True
    },
    {
        'file': r'G:\My Drive\Research\PVRD1\Literature\PID_degradation_time\Oh-MicroelectronicsReliability_2017-Fig1_Pmax_1000V_85C_85RH.csv',
        'time_units': 'min',
        'label': 'Oh 2017',
        'color': 'tab:orange',
        'marker': 's',
        'type': 'power',
        'normalized': True
    }
]

csv_index = r'file_index.csv'


defaultPlotStyle = {'font.size': 11,
                    'font.family': 'Arial',
                    'font.weight': 'regular',
                    'legend.fontsize': 12,
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
                    'xtick.top': False,
                    'ytick.right': False,
                    'lines.linewidth': 2.5,
                    'lines.markersize': 10,
                    'lines.markeredgewidth': 0.85,
                    'axes.labelpad': 5.0,
                    'axes.labelsize': 12,
                    'axes.labelweight': 'regular',
                    'legend.handletextpad': 0.2,
                    'legend.borderaxespad': 0.2,
                    'axes.linewidth': 1.25,
                    'axes.titlesize': 12,
                    'axes.titleweight': 'bold',
                    'axes.titlepad': 6,
                    'figure.titleweight': 'bold',
                    'figure.dpi': 100}

if __name__ == '__main__':
    mpl.rcParams.update(defaultPlotStyle)
    xfmt = ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((-3, 3))

    
        
    if platform.system() == 'Windows':
        base_folder = r'\\?\\' + os.path.abspath(base_folder)
        results_folder = r'\\?\\' + results_folder

    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
        
    n_plots = len(simulation_result_dirs)
    failure_times = np.empty(
        n_plots, dtype=np.dtype([
            (sweep_variable, 'd'),
            ('t 5% loss (h)', 'd'),
            ('t 10% loss (h)', 'd'),
            ('t 15% loss (h)', 'd'),
            ('t 20% loss (h)', 'd'),
        ])
    )
    
    fig = plt.figure()
    fig.set_size_inches(5.0, 6.0, forward=True)
    fig.subplots_adjust(hspace=0.2, wspace=0.35)
    gs0 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, width_ratios=[1])
    gs00 = gridspec.GridSpecFromSubplotSpec(nrows=3, ncols=1, subplot_spec=gs0[0])
    ax1 = fig.add_subplot(gs00[0, 0])
    ax2 = fig.add_subplot(gs00[1, 0])
    ax3 = fig.add_subplot(gs00[2, 0])
    
    
    c_map1 = mpl.cm.get_cmap('winter')
    c_map2 = mpl.cm.get_cmap('spring')
    c_map3 = mpl.cm.get_cmap('summer')
    c_map4 = mpl.cm.get_cmap('autumn')
    
    if sweep_log:
        normalize = mpl.colors.LogNorm(vmin=min(sweep), vmax=max(sweep))
    else:
        normalize = mpl.colors.Normalize(vmin=0.0, vmax=t_max)
    colors_1 = [c_map1(normalize(x)) for x in sweep]
    scalar_maps_1 = mpl.cm.ScalarMappable(cmap=c_map1, norm=normalize)
    colors_2 = [c_map2(normalize(x)) for x in sweep]
    scalar_maps_2 = mpl.cm.ScalarMappable(cmap=c_map2, norm=normalize)
    colors_3 = [c_map3(normalize(x)) for x in sweep]
    scalar_maps_3 = mpl.cm.ScalarMappable(cmap=c_map3, norm=normalize)
    colors_4 = [c_map4(normalize(x)) for x in sweep]
    scalar_maps_4 = mpl.cm.ScalarMappable(cmap=c_map4, norm=normalize)
        
    
    for i, fn, sv in zip(range(n_plots), simulation_result_dirs, sweep):
        eff_df = pd.read_csv(
            filepath_or_buffer=os.path.join(
                base_folder, fn, 'jv_plots', 'efficiency_results.csv'
            )
        )
        rsh_df = pd.read_csv(
            filepath_or_buffer=os.path.join(
                base_folder, fn, 'analysis_plots', 'rsh_data.csv'
            )
        )
        
        time = eff_df['time (s)'].to_numpy() / 3600
        jsc = eff_df['jsc (mA/cm2)'].to_numpy()
        voc = eff_df['voc (V)'].to_numpy()
        pmpp = eff_df['pd_mpp (mW/cm2)'].to_numpy()
        ff = pmpp/jsc/voc
        rsh = rsh_df['Rsh (Ohms cm2)'].to_numpy()
        
        L_diff = np.sqrt(sv*eff_df['time (s)'].to_numpy())
        
        
        ax1.plot(
            time, pmpp/pmpp[0], color=c_map1(normalize(sv)), ls='-',
            # marker='o', fillstyle='none'
        )
        ax2.plot(
            time, rsh, color=c_map1(normalize(sv)), ls='-',
            # marker='o', fillstyle='none'
        )
        ax3.plot(
            time, ff, color=c_map1(normalize(sv)), ls='-',
            # marker='o', fillstyle='none'
        )
        
    for i, lf in enumerate(literature_files):
        fn = lf['file']
        time_units = lf['time_units']
        label = lf['label']
        if lf['type'] == 'power':
            column_names = ['time', 'power']
        else:
            column_names = ['time', 'Rsh']
        lit_df = pd.read_csv(fn, skiprows=1, header=0, names=column_names)
        
        if time_units == 'min':
            time_lf = lit_df['time'].to_numpy()/60
        elif time_units == 's':
            time_lf = lit_df['time'].to_numpy()/3600
        elif time_units == 'h':
            time_lf = lit_df['time'].to_numpy()
        else:  # Assume hours
            time_lf = lit_df['time'].to_numpy()
        
        
        if lf['type'] == 'power':
            normalized_power = lit_df['power']
            ax1.plot(
                time_lf, normalized_power, 'o-', fillstyle='none', 
                color=lf['color'], label=label,
                marker=lf['marker']
            )
        else:
            rsh = lit_df['Rsh']
            if not lf['normalized']:
                rsh = rsh/rsh[0]
            ax2.plot(
                time_lf, rsh, 's-', fillstyle='none', 
                color=lf['color'], label=label,
                marker=lf['marker']
            )
        
    
    ax1.set_ylabel('$P_{\\mathrm{mpp}}/P_{\\mathrm{mpp}}(t=0)$')
    ax2.set_ylabel('$R_{\mathrm{sh}}$ ($\Omega$ cm$^{2}$)')
    ax3.set_ylabel('FF')
    ax3.set_xlabel('Time (hr)')
    
    ax1.set_xlim(0, t_max)
    ax2.set_xlim(0, t_max)
    ax3.set_xlim(0, t_max)
    
    ax1.tick_params(labelbottom=False, top=True, right=True, which='both', labeltop=True)
    ax2.tick_params(labelbottom=False, top=False, right=True, which='both')
    ax3.tick_params(labelbottom=True, top=False, right=True, which='both')
    
    
    ax2.set_yscale('log')
    ax2.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=5))
    ax2.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=50, subs=np.arange(2, 10) * .1))

    ax1.xaxis.set_major_formatter(xfmt)
    ax1.xaxis.set_major_locator(mticker.MaxNLocator(12, prune=None))
    ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))

    ax1.yaxis.set_major_formatter(xfmt)
    ax1.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    ax1.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))

    ax2.xaxis.set_major_formatter(xfmt)
    ax2.xaxis.set_major_locator(mticker.MaxNLocator(12, prune=None))
    ax2.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))


    ax3.xaxis.set_major_formatter(xfmt)
    ax3.xaxis.set_major_locator(mticker.MaxNLocator(12, prune=None))
    ax3.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))

    ax3.yaxis.set_major_formatter(xfmt)
    ax3.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    ax3.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="5%", pad=0.03)
    cbar1 = fig.colorbar(scalar_maps_1, cax=cax1)
    cbar1.set_label(sweep_variable, rotation=90)
    
    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="5%", pad=0.03)
    cbar2 = fig.colorbar(scalar_maps_1, cax=cax2)
    cbar2.set_label(sweep_variable, rotation=90)
    
    divider3 = make_axes_locatable(ax3)
    cax3 = divider3.append_axes("right", size="5%", pad=0.03)
    cbar3 = fig.colorbar(scalar_maps_1, cax=cax3)
    cbar3.set_label(sweep_variable, rotation=90)
    
    lg1 = ax1.legend(frameon=False, fontsize=12, loc='upper right')
    # lg2 = ax2.legend(frameon=False, fontsize=12, loc='upper right')

    plt.tight_layout()
    fig.savefig(os.path.join(results_folder, 'time_series.png'), dpi=600)