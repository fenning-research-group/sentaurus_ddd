# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 07:56:11 2020

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
t_max = 12
simulation_result_dirs = [
        r'Cs=1E+16D1=4E-16_h=1E-12_D2=1E-14_rho=4E-05_SD=1.0_s=5E+01',
        r'Cs=1E+16D1=4E-16_h=1E-08_D2=1E-14_rho=4E-05_SD=1.0_s=5E+01',
        # r'Cs=1E+16D1=4E-16_h=1E-04_D2=1E-14_rho=4E-05_SD=1.0_s=5E+01',
    ]

sweep_variable = r'$h$ (cm/s)'
sweep_variable_units = r'(cm/s)'
swee_variable_name = r'$h$'
sweep = [1E-12, 1E-8]
sweep_log = True
label_literature = True

results_folder = r'G:\Shared drives\FenningLab2\Projects\PVRD1\Simulations\Sentaurus PID\results\analysis\0.5MVcm_D2=1E-14_Cs=1E16_variable_h'

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
        'file': r'G:\My Drive\Research\PVRD1\Literature\PID_degradation_time\Hacke_ProgressInPhoto2013_Fig1_600V_85C_type2_1.csv',
        'time_units': 'h',
        'label': 'Hacke 2013',
        'color': 'tab:orange',
        'marker': 's',
        'type': 'power',
        'normalized': True
    },
    {
        'file': r'G:\My Drive\Research\PVRD1\Literature\PID_degradation_time\Lausch_IEEEJPV_2014_600V_Rsh.csv',
        'time_units': 'h',
        'label': 'Lausch 2014',
        'color': 'tab:red',
        'marker': '^',
        'type': 'Rsh',
        'normalized': False
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
    fig.set_size_inches(5.0, 5.0, forward=True)
    fig.subplots_adjust(hspace=0.2, wspace=0.35)
    gs0 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, width_ratios=[1])
    gs00 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs0[0])
    ax1 = fig.add_subplot(gs00[0, 0])
    ax2 = fig.add_subplot(gs00[1, 0])
    # ax3 = fig.add_subplot(gs00[2, 0])
    
    
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
    
    zorder = 0
    for i, lf in enumerate(literature_files):
        fn = lf['file']
        time_units = lf['time_units']
        label = lf['label'] if label_literature else None
        if lf['type'] == 'power':
            column_names = ['time', 'power']
        else:
            column_names = ['time', 'Rsh']
        lit_df = pd.read_csv(fn, skiprows=0, header=0, names=column_names)
        
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
                time_lf, normalized_power, fillstyle='none', ls='--',
                color=lf['color'], label=label, marker=lf['marker'],
                zorder=zorder
            )
        else:
            rsh = lit_df['Rsh']
            if not lf['normalized']:
                rsh = rsh/rsh[0]
            ax2.plot(
                time_lf, rsh, fillstyle='none', ls='--',
                color=lf['color'], label=label, marker=lf['marker'],
                zorder=zorder
            )
        zorder += 1
        
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
        
        t_interp = np.linspace(np.amin(time), np.amax(time), num=200)
        f_p_interp = interpolate.interp1d(time, pmpp, kind='linear')
        pmpp_interp = f_p_interp(t_interp)
        
        idx_5 = (np.abs(pmpp_interp/pmpp_interp[0] - 0.95)).argmin()
        idx_10 = (np.abs(pmpp_interp/pmpp_interp[0] - 0.9)).argmin()
        idx_15 = (np.abs(pmpp_interp/pmpp_interp[0] - 0.85)).argmin()
        idx_20 = (np.abs(pmpp_interp/pmpp_interp[0] - 0.8)).argmin()
        
        failure_times[i] = (
            sv,
            t_interp[idx_5], t_interp[idx_10], t_interp[idx_15], 
            t_interp[idx_20], 
        )
        
        L_diff = np.sqrt(sv*eff_df['time (s)'].to_numpy())
        
        sv_str = '{0:.1E}'.format(sv)
        sv_arr = sv_str.split('E')
        sv_arr = np.array(sv_arr,dtype=float)
        sv_txt = r'{0} = $\mathregular{{ 10^{{{1:.0f}}} }}$ {2}'.format(swee_variable_name, sv_arr[1], sweep_variable_units)
        ax1.plot(
            time, pmpp/pmpp[0], color=c_map1(normalize(sv)), ls='-',
            zorder=zorder, label=sv_txt
            # marker='o', fillstyle='none'
        )
        ax2.plot(
            time, rsh/rsh[0], color=c_map1(normalize(sv)), ls='-',
            zorder=zorder, #label=sv_txt
            # marker='o', fillstyle='none'
        )
        # ax3.plot(
        #     time, ff, color=c_map1(normalize(sv)), ls='-',
        #     zorder=zorder
        #     # marker='o', fillstyle='none'
        # )
        zorder += 1
    
        
    
    ax1.set_ylabel('$P_{\\mathrm{mpp}}/P_{\\mathrm{mpp}}(t=0)$')
    # ax2.set_ylabel('$R_{\mathrm{sh}}$ ($\Omega$ cm$^{2}$)')
    ax2.set_ylabel('$R_{\mathrm{sh}}/R_{\mathrm{sh}}(t=0)$')
    # ax3.set_ylabel('FF')
    ax2.set_xlabel('Time (hr)')
    
    ax1.set_xlim(0, t_max)
    ax2.set_xlim(0, t_max)
    # ax3.set_xlim(0, t_max)
    
    ax1.tick_params(labelbottom=False, top=True, right=True, which='both', labeltop=True)
    # ax2.tick_params(labelbottom=False, top=False, right=True, which='both')
    ax2.tick_params(labelbottom=True, top=False, right=True, which='both')
    
    
    ax2.set_yscale('log')
    ax2.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=5))
    ax2.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=50, subs=np.arange(2, 10) * .1))

    ax1.xaxis.set_major_formatter(xfmt)
    ax1.xaxis.set_major_locator(mticker.MaxNLocator(6, prune=None))
    ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))

    ax1.yaxis.set_major_formatter(xfmt)
    ax1.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    ax1.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))

    ax2.xaxis.set_major_formatter(xfmt)
    ax2.xaxis.set_major_locator(mticker.MaxNLocator(6, prune=None))
    ax2.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))


    # ax3.xaxis.set_major_formatter(xfmt)
    # ax3.xaxis.set_major_locator(mticker.MaxNLocator(12, prune=None))
    # ax3.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))

    # ax3.yaxis.set_major_formatter(xfmt)
    # ax3.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    # ax3.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    # divider1 = make_axes_locatable(ax1)
    # cax1 = divider1.append_axes("right", size="5%", pad=0.03)
    # cbar1 = fig.colorbar(scalar_maps_1, cax=cax1)
    # cbar1.set_label(sweep_variable, rotation=90)
    
    # divider2 = make_axes_locatable(ax2)
    # cax2 = divider2.append_axes("right", size="5%", pad=0.03)
    # cbar2 = fig.colorbar(scalar_maps_1, cax=cax2)
    # cbar2.set_label(sweep_variable, rotation=90)
    
    # divider3 = make_axes_locatable(ax3)
    # cax3 = divider3.append_axes("right", size="5%", pad=0.03)
    # cbar3 = fig.colorbar(scalar_maps_1, cax=cax3)
    # cbar3.set_label(sweep_variable, rotation=90)
    
    lg1 = ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    if label_literature:
        lg2 = ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)


    plt.tight_layout()
    fig.savefig(os.path.join(results_folder, 'time_series.png'), dpi=600)
    
    
    """
    Plot loss as a function of electric field
    """
    
    fig = plt.figure()
    fig.set_size_inches(4.0, 3.5, forward=True)
    # fig.subplots_adjust(hspace=0.2, wspace=0.35)
    gs0 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, width_ratios=[1])
    gs00 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs0[0])
    ax1 = fig.add_subplot(gs00[0, 0])
    
    idx5_neq_0 = failure_times['t 5% loss (h)'] > 0
    idx10_neq_0 = failure_times['t 10% loss (h)'] > 0
    idx15_neq_0 = failure_times['t 15% loss (h)'] > 0
    idx20_neq_0 = failure_times['t 20% loss (h)'] > 0
    
    ax1.set_xscale('log')
    # ax1.set_yscale('log')
    
    ax1.plot(
        failure_times[sweep_variable][idx5_neq_0], failure_times['t 5% loss (h)'][idx5_neq_0],
        ls='--', marker='o', label='5% loss', fillstyle='none'
    )
    
    ax1.plot(
        failure_times[sweep_variable][idx10_neq_0], failure_times['t 10% loss (h)'][idx10_neq_0],
        ls='--', marker='o', label='10% loss', fillstyle='none'
    )
    
    ax1.plot(
        failure_times[sweep_variable][idx15_neq_0], failure_times['t 15% loss (h)'][idx15_neq_0],
        ls='--', marker='o', label='15% loss', fillstyle='none'
    )
    
    ax1.plot(
        failure_times[sweep_variable][idx20_neq_0], failure_times['t 20% loss (h)'][idx20_neq_0],
        ls='--', marker='o', label='20% loss', fillstyle='none'
    )
    
    ax1.set_xlabel(sweep_variable)
    ax1.set_ylabel('Failure Time (h)')
    ax1.set_ylim(top=3)
    
    ax1.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=5))
    ax1.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=50, subs=np.arange(2, 10) * .1))
    
    ax1.legend(
        # bbox_to_anchor=(0., 1.02, 1., .102), 
        # loc='lower left', ncol=2, mode='expand', borderaxespad=0.
        loc='upper right', ncol=2, mode=None, borderaxespad=0.,
        frameon=False
    )
    
    ax1.set_title('$\mathregular{C_s = 10^{16}\;cm^{-3}}$, $\mathregular{D_{Si}}$ = $\mathregular{10^{-14}\;cm^2/s}$, 85 °C')
    plt.tight_layout()
    fig.savefig(os.path.join(results_folder, 'failure_time_Cs=1E16_0.5MVcm_varying_h_85C.png'), dpi=600)
    
    df_degradation = pd.DataFrame(failure_times)
    df_degradation.to_csv(
        path_or_buf=os.path.join(
            results_folder,
            'failure_time_Cs=1E16_0.5MVcm_varying_h_85C.csv'
        ),
        index=False
    )