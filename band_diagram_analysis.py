# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 06:47:30 2020

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

base_folder = r'C:\Users\Erick\PycharmProjects\sentaurus_ddd\results\work_function'
output_folder = r'all_plots'
csv_index = r'C:\Users\Erick\PycharmProjects\sentaurus_ddd\results\work_function\all_plots\band_diagrams.csv'
shunt_depth = 1.0

xfmt = ScalarFormatter(useMathText=True)
xfmt.set_powerlimits((-3, 3))

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


cmap_blues = mpl.cm.get_cmap('Blues')
cmap_greens = mpl.cm.get_cmap('Greens')
cmap_greys = mpl.cm.get_cmap('Greys')
cmap_purples = mpl.cm.get_cmap('Purples')
cmap_reds = mpl.cm.get_cmap('Reds')


if __name__ == '__main__':
    mpl.rcParams.update(defaultPlotStyle)
    if platform.system() == 'Windows':
        base_folder = r'\\?\\' + os.path.abspath(base_folder)
        csv_index = r'\\?\\' + csv_index
    # Read the csv index
#    shunt_depths = index_df['shunt depth (um)'].unique()
    index_df = pd.read_csv(filepath_or_buffer=csv_index)
    r0 = index_df[index_df['shunt depth (um)'] == 0.0]
    fn0 = r0['filename'][0]
    work_function0 = float(r0['work_function'])
    bd0 = pd.read_csv(os.path.join(base_folder, fn0), skiprows=[1,2])
    bd0 = bd0[bd0['Y'] <= 1.0].reset_index(drop=True)
    bd0 = bd0[bd0['Y'] >=0 ].reset_index(drop=True)
    
    
    index_df = index_df[index_df['shunt depth (um)'] == shunt_depth].reset_index(drop=True)

    
    
    work_functions = np.array(index_df['work_function'])
    normalize = mpl.colors.Normalize(vmin=np.amin(work_functions), vmax=np.amax(work_functions))
    colors = [cmap_greens(normalize(c)) for c in work_functions]
    scalar_maps = mpl.cm.ScalarMappable(cmap=cmap_greens, norm=normalize)
    
    

    #plt.ioff()
    fig = plt.figure()
    fig.set_size_inches(6.0, 3.0, forward=True)
    fig.subplots_adjust(hspace=0.35, wspace=0.55)
    gs0 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, width_ratios=[1])
    gs00 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs0[0])
    ax1 = fig.add_subplot(gs00[0, 0])
    ax2 = fig.add_subplot(gs00[0, 1])

    for j, r in index_df.iterrows():
        fn = r['filename']
        wf = float(r['work_function'])
        print('Analyzing file: {0}'.format(fn))
        bd = pd.read_csv(os.path.join(base_folder, fn), skiprows=[1,2])
        bd = bd[bd['Y'] <= 1.75].reset_index(drop=True)
        bd = bd[bd['Y'] >=0.75].reset_index(drop=True)
        # Plot Conduction Band
        ax2.plot(bd['Y'], bd['ConductionBandEnergy'], color=cmap_reds(normalize(wf)))
        # Plot Valence Band
        ax2.plot(bd['Y'], bd['ValenceBandEnergy'], color=cmap_blues(normalize(wf)))
        # Plot EFn
        ax2.plot(bd['Y'], bd['eQuasiFermiEnergy'], '--' , dashes=(3, 1),  
                 color=cmap_purples(normalize(wf)))
        # Plot EFp
        ax2.plot(bd['Y'], bd['hQuasiFermiEnergy'], '--', dashes=(3, 1), 
                 color=cmap_greens(normalize(wf)))
        
    # Plot Conduction Band
    ax1.plot(bd0['Y'], bd0['ConductionBandEnergy'], color=cmap_reds(normalize(wf)),
             label='$E_c$')
    # Plot Valence Band
    ax1.plot(bd0['Y'], bd0['ValenceBandEnergy'], color=cmap_blues(normalize(wf)),
             label='$E_v$')
    # Plot EFn
    ax1.plot(bd0['Y'], bd0['eQuasiFermiEnergy'], '--' , dashes=(3, 1), 
             color=cmap_purples(normalize(wf)),
             label='$E_{Fn}$')
    # Plot EFp
    ax1.plot(bd0['Y'], bd0['hQuasiFermiEnergy'], '--' , dashes=(3, 1),
             color=cmap_greens(normalize(wf)),
             label='$E_{Fp}$')


    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.03)
    cbar = fig.colorbar(scalar_maps, cax=cax)
    cbar.set_label('Work Function (eV)', rotation=90)

    ax1.set_xlabel('Depth (um)')
    ax1.set_ylabel('Energy (eV)')
    
    ax2.set_xlabel('Depth (um)')
    ax2.set_ylabel('Energy (eV)')
    
    ax1.legend(loc='upper left', frameon=False, fontsize=11)
    
    
    ax1.xaxis.set_major_formatter(xfmt)
    ax1.xaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    ax2.xaxis.set_major_formatter(xfmt)
    ax2.xaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    ax2.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))

    ax1.yaxis.set_major_formatter(xfmt)
    ax1.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    ax1.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    ax2.yaxis.set_major_formatter(xfmt)
    ax2.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    ax2.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    ax1.set_ylim(top=1.25)
    
    

    ax1.set_title('Un-shunted')
    ax2.set_title('Shunt Depth: {0:.1f} um'.format(shunt_depth))
    
    

    plt.tight_layout()
    # plt.show()
    fig.savefig(os.path.join(base_folder, output_folder, 'band_diagram_sd_{0:.1f}_um.png'.format(shunt_depth)), dpi=600)



