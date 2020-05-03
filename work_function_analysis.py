import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import ScalarFormatter
import platform
import os
import pid.analysis as pia
import matplotlib.gridspec as gridspec
from scipy import interpolate

base_folder = r'C:\Users\Erick\PycharmProjects\sentaurus_ddd\results\work_function'
output_folder = r'all_plots'
csv_index = r'C:\Users\Erick\PycharmProjects\sentaurus_ddd\results\work_function\all_plots\file_index.csv'
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

nu_t_dtype = np.dtype([
    ('shunt depth (um)', 'd'),
    ('conductivity (Ohm/cm)^{-1}', 'd'),
    ('jsc (mA/cm2)', 'd'),
    ('voc (V)', 'd'),
    ('pd_mpp (mW/cm2)', 'd'),
    ('v_mpp (V)', 'd'),
    ('j_mpp (mA/cm2)', 'd'),
    ('efficiency', 'd')
])

if __name__ == '__main__':
    mpl.rcParams.update(defaultPlotStyle)
    if platform.system() == 'Windows':
        base_folder = r'\\?\\' + os.path.abspath(base_folder)
        csv_index = r'\\?\\' + csv_index
    # Read the csv index
    index_df = pd.read_csv(filepath_or_buffer=csv_index)
    shunt_depths = index_df['shunt depth (um)'].unique()
    # A color palette for each depth
    color_palettes = ['Blues', 'Oranges']

    r0 = index_df.iloc[0]
    fn0 = r0['filename']
    work_function0 = float(r0['work_function'])
    analyzer = pia.Analysis(folder_path=base_folder)
    jv_curve0 = analyzer.read_jv(h5_filename=os.path.join(base_folder, fn0))
    nu_data0 = analyzer.estimate_efficiency(voltage=jv_curve0['voltage (V)'],
                                            current=jv_curve0['current (mA/cm2)'])
    
    jsc_max = nu_data0['jsc']

    idx_curve = jv_curve0['current (mA/cm2)'] >= 0
    voltage0 = jv_curve0['voltage (V)'][idx_curve]
    current0 = jv_curve0['current (mA/cm2)'][idx_curve]
    efficiency0 = analyzer.estimate_efficiency(voltage0, current0)

    for i, sd in enumerate(shunt_depths[shunt_depths > 0]):
        cmap = mpl.cm.get_cmap(color_palettes[i])
        files_df: pd.DataFrame = index_df[index_df['shunt depth (um)'] == sd].reset_index(drop=True)
        work_functions = np.array(files_df['work_function'])
        normalize = mpl.colors.Normalize(vmin=np.amin(work_functions), vmax=np.amax(work_functions))
        colors = [cmap(normalize(c)) for c in work_functions]
        scalar_maps = mpl.cm.ScalarMappable(cmap=cmap, norm=normalize)
        results = np.empty(len(files_df), dtype=nu_t_dtype)
        

        #plt.ioff()
        fig = plt.figure()
        fig.set_size_inches(6.0, 3.0, forward=True)
        fig.subplots_adjust(hspace=0.35, wspace=0.65)
        gs0 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, width_ratios=[1])
        gs00 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs0[0])
        ax1 = fig.add_subplot(gs00[0, 0])
        ax2 = fig.add_subplot(gs00[0, 1])

        for j, r in files_df.iterrows():
            fn = r['filename']
            wf = float(r['work_function'])
            print('Analyzing file: {0}'.format(fn))
            jv_curve = analyzer.read_jv(h5_filename=os.path.join(base_folder, fn))
            nu_data = analyzer.estimate_efficiency(voltage=jv_curve['voltage (V)'],
                                                   current=jv_curve['current (mA/cm2)'])
            results[j] = (
                sd, wf, nu_data['jsc'], nu_data['voc'], nu_data['pd_mpp'], nu_data['v_mpp'], nu_data['j_mpp'],
                nu_data['efficiency']
            )
            idx_curve = jv_curve['current (mA/cm2)'] >= 0
            voltage = jv_curve['voltage (V)'][idx_curve]
            current = jv_curve['current (mA/cm2)'][idx_curve]
            ax1.plot(voltage, current, color=cmap(normalize(wf)))

        ax1.plot(voltage0, current0, color='r')

        # interpoate the efficiency
        # interpolation = interpolate.splrep(conductivities, results['efficiency'], k=3)
        f = interpolate.interp1d(work_functions, results['efficiency'], kind='slinear', fill_value="extrapolate")
        t_interp = np.linspace(np.amin(wf), np.amax(wf), num=1000)
        nu_interp = f(t_interp)
        # nu_interp = interpolate.splev(t_interp, interpolation)

        ax2.plot(work_functions, results['efficiency'] * 100, marker='o', ls='none', fillstyle='none',
                 color=cmap(normalize(np.amax(wf))),
                 label='Simulation')

        ax2.plot(work_function0, efficiency0['efficiency']*100, marker='o', color='r', ls='none', fillstyle='none')

        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.03)
        cbar = fig.colorbar(scalar_maps, cax=cax)
        cbar.set_label('Work Function (S/cm)', rotation=90)

        ax1.set_xlabel('Bias (V)')
        ax1.set_ylabel('J (mA/cm$^2$)')
        
#        ax1.set_ylabel('J$_{\mathregular{sc}}$ unshunted - J (mA/cm$^2$)')

        ax2.set_xlabel('Work Function (S/cm)')
        ax2.set_ylabel('Efficiency (%)')

#        ax2.plot(t_interp, nu_interp * 100, ':', color='k', dashes=(3, 1), label='Guide to the eye', lw=1.5)
#        ax1.set_yscale('log')
#        ax1.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=6))
#        ax1.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=60, subs=np.arange(2, 10) * .1))

        ax2.yaxis.set_major_formatter(xfmt)
        ax2.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
        ax2.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
#        ax2.set_xscale('log')
#        ax2.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=6))
#        ax2.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=60, subs=np.arange(2, 10) * .1))

        ax1.set_title('Shunt Depth: {0:.1f} um'.format(sd))
        ax2.set_title('Shunt Depth: {0:.1f} um'.format(sd))
        
        ax1.text(0.05, 0.15, 'Unshunted',
                 horizontalalignment='left',
                 verticalalignment='bottom', 
                 transform=ax1.transAxes,
                 fontsize=12,
                 color='r')
        
        ax1.text(0.05, 0.05, 'Shunted',
                 horizontalalignment='left',
                 verticalalignment='bottom', 
                 transform=ax1.transAxes,
                 fontsize=12,
                 color='tab:blue')

        plt.tight_layout()
        # plt.show()
        fig.savefig(os.path.join(base_folder, output_folder, 'efficiency_sd_{0:.1f}_um.png'.format(sd)), dpi=600)

