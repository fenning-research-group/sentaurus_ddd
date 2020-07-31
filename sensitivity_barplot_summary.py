# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 06:21:36 2020

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

results_folder = r'G:\Shared drives\FenningLab2\Projects\PVRD1\Simulations\Sentaurus PID\results\analysis\sensitivity_barplot_summary' 

data_files = [
    {
        'filepath': r'G:\Shared drives\FenningLab2\Projects\PVRD1\Simulations\Sentaurus PID\results\analysis\Cs1E18_E=0.01_variable_D2\failure_time_10kVcm_85C.csv',
        'sweep_variable' : '$D_{\mathrm{Si}}$',
        'sweep_variable_units' : '$\mathrm{cm^{-3}}$',
        'unit_prefix': None,
        'columns' : ['D (cm^2/s)', '5% loss (h)'],
        'sweep_values' : [1e-16, 1e-14],
    },
    {
         'filepath' : r'G:\Shared drives\FenningLab2\Projects\PVRD1\Simulations\Sentaurus PID\results\analysis\0.5MVcm_D2=1E-14_h=1E-12_variable_Cs\failure_times_D2=1E-14_h=1E-12_different_cs.csv',
         'sweep_variable' : '$C_S$',
         'sweep_variable_units' : '$\mathrm{cm^{-3}}$',
         'columns' : ['$C_s$ (cm$\mathregular{^{-3}}$)', '5% loss (h)'],
         'sweep_values' : [1e+16, 1e+22],
    },
    {
         'filepath' : r'G:\Shared drives\FenningLab2\Projects\PVRD1\Simulations\Sentaurus PID\results\analysis\Cs1E16_D2=1E-14_h=1E-12_variable_E\failure_times_D2=1E-14_h=1E-12_different_cs.csv',
         'sweep_variable' : '$E$',
         'sweep_variable_units' : 'V/cm',
         'unit_prefix' : 'M',
         'columns' : ['$E$ (MV/cm)', '5% loss (h)'],
         'sweep_values' : [0.01, 0.5],
    },
    {
         'filepath' : r'G:\Shared drives\FenningLab2\Projects\PVRD1\Simulations\Sentaurus PID\results\analysis\0.5MVcm_D2=1E-14_Cs=1E16_variable_h\failure_time_Cs=1E16_0.5MVcm_varying_h_85C.csv',
         'sweep_variable' : '$h$',
         'sweep_variable_units' : 'cm/s',
         'unit_prefix' : None,
         'columns' : ['$h$ (cm/s)', '5% loss (h)'],
         'sweep_values' : [1E-12, 1E-8],
    }
    
]

bar_width = 0.25

def autolabel(rects, data_files_, column_index, xpos='center'):
    """
    Attach a text label above each bar in *rects*, displaying its height.

    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """

    xpos = xpos.lower()  # normalize the case of the parameter
    ha = {'center': 'center', 'right': 'left', 'left': 'right'}
    offset = {'center': 0.5, 'right': 0.57, 'left': 0.43}  # x_txt = x + w*off
        
    
    for i, rect in enumerate(rects):
        height = rect.get_height()
        fi = data_files_[i]
        
        sv = fi['sweep_values'][column_index]
        val_str = '{0:.1E}'.format(sv)
        val_arr = val_str.split('E')
        val_arr = np.array(val_arr,dtype=float)
        lbl = r' $\mathregular{{ 10^{{{0:.0f}}} }}$ {1}'.format(
            val_arr[1], fi['sweep_variable_units']
        )
        if fi['sweep_variable'] == '$E$':
            if sv * 10 > 1:
                lbl = ' {0:.1f} MV/cm'.format(sv)
            elif sv < 0.1:
                lbl = ' {0:.0f} kV/cm'.format(sv*1000)
        if height < 2:
            ax1.text(
                rect.get_x() + rect.get_width()*offset[xpos], 1.01*height, lbl, 
                ha=ha[xpos], va='bottom', rotation=90, fontsize=9
            )
        else:
            ax1b.text(
                rect.get_x() + rect.get_width()*offset[xpos], 1.01*height,lbl, 
                ha=ha[xpos], va='bottom', rotation=90, fontsize=9 
            )

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
        results_folder = r'\\?\\' + results_folder

    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    
    fig = plt.figure()
    fig.set_size_inches(5.0,3.75, forward=True)
    fig.subplots_adjust(hspace=0.1, wspace=0.4)
    gs0 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, width_ratios=[1])
    gs00 = gridspec.GridSpecFromSubplotSpec(
        nrows=2, ncols=1, height_ratios=[1, 1],
        subplot_spec = gs0[0]
    )
    n_variables = len(data_files)
    sweep_values = np.empty(
        n_variables, dtype=np.dtype([('low','d'), ('high','d')])
    )
    xticks = list(range(n_variables))
    xtick_labels = []
    ax1 = fig.add_subplot(gs00[1,0])
    ax1b = fig.add_subplot(gs00[0,0], sharex=ax1)
    
    ax1.set_ylim(top=2.5)
    ax1b.set_ylim(bottom=30, top=45)
    
    
    zorder = 0
    for i, fi in enumerate(data_files):
        fn = fi['filepath']
        if platform.system() == 'Windows':
            fn = r'\\?\\' + fn
        xtick_labels.append(
            '{0}'.format(fi['sweep_variable'])
        )
        df = pd.read_csv(filepath_or_buffer=fn, usecols=[0,1] )
        df_cols = df.columns
        values = df[df.iloc[:,0].isin(fi['sweep_values'])].sort_values(by=[df_cols[0]])
        sweep_values[i] =tuple(values.iloc[:,1])
        
        print(values)
        
    rects1 = ax1.bar(
        xticks, sweep_values['low'], color='tab:red', label='Low', 
        width=bar_width,
    )
    rects2 = ax1.bar(
        [x +bar_width for x in xticks], sweep_values['high'], 
        color='tab:green', label='High', width=bar_width,
    )
    
    rects1b = ax1b.bar(
        xticks, sweep_values['low'], color='tab:red', label='Low', 
        width=bar_width,
    )
    rects2b = ax1b.bar(
        [x +bar_width for x in xticks], sweep_values['high'], 
        color='tab:green', label='High', width=bar_width,
    )
    
    
        
    
    ax1.spines['top'].set_visible(False)
    ax1b.spines['bottom'].set_visible(False)
    
    ax1.xaxis.tick_bottom()
    ax1b.xaxis.tick_top()
    ax1b.tick_params(labeltop=False)
    
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1b.transAxes, color='k', clip_on=False)
    ax1b.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax1b.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    
    kwargs.update(transform=ax1.transAxes)  # switch to the bottom axes
    ax1.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xtick_labels)
    
    ax1.yaxis.set_major_formatter(xfmt)
    ax1.yaxis.set_major_locator(mticker.MaxNLocator(5, prune='upper'))
    ax1.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    ax1b.yaxis.set_major_formatter(xfmt)
    ax1b.yaxis.set_major_locator(mticker.MaxNLocator(4, prune='lower'))
    ax1b.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    fig.text(
        0.0, 0.5, '5% Failure Time (h)', va='center', rotation='vertical',
        fontsize=12
    )
    
    ax1.set_xlabel('Simulation parameter')
    
    autolabel(rects1, data_files, 0, "center")
    autolabel(rects2, data_files, 1, "center")
    fig.savefig(os.path.join(results_folder, 'failure_time_85C_summary.png'), dpi=600)