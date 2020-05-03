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
import pid.confidence as cf
import matplotlib.gridspec as gridspec
from scipy import interpolate
import h5py
from scipy import optimize
from scipy.linalg import svd

# radius = 50E-4 # cm
area = np.pi*(5)*1E-6


base_folder = r'C:\Users\Erick\PycharmProjects\sentaurus_ddd\results\time_dependence\20200428_sc=1E2-D2=1E-14_shuntDepth=1.0-0'
output_folder = r'analysis_plots'
csv_index = r'file_index.csv'
na_file = r'C:\Users\Erick\PycharmProjects\sentaurus_ddd\results\time_dependence\20200428_sc=1E2-D2=1E-14_shuntDepth=1.0-0\two_layers_D1=4E-16cm2ps_D2=1E-14cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-10_m1.0e+00_pnp.h5'
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
    ('time (s)', 'd'),
    ('jsc (mA/cm2)', 'd'),
    ('voc (V)', 'd'),
    ('pd_mpp (mW/cm2)', 'd'),
    ('v_mpp (V)', 'd'),
    ('j_mpp (mA/cm2)', 'd'),
    ('efficiency', 'd')
])

rsh_dtype = np.dtype([
            ('time (h)', 'd'),
            ('Rsh (Ohms cm2)', 'd'),
            ('Rsh lci', 'd'),
            ('Rsh uci', 'd')
        ])

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def fobj(a: np.ndarray, X: np.ndarray, Y: np.ndarray):
#    a = np.power(10, a)
    return -a[0]*X + a[1] - Y

def jac(a: np.ndarray, X: np.ndarray, Y: np.ndarray):
#    a = np.power(10, a)
    n = len(X)
    res = np.empty((n, 2))
    for i, x in enumerate(X):
        res[i] = (-x, 1.0)
    return res

def find_pcov(optimize_result):
    popt = optimize_result.x
    ysize = len(optimize_result.fun)
    cost = 2 * optimize_result.cost  # res.cost is half sum of squares!
    s_sq = cost / (ysize - popt.size)
    
    # Do Moore-Penrose inverse discarding zero singular values.
    _, s, VT = svd(optimize_result.jac, full_matrices=False)
    threshold = np.finfo(float).eps * max(optimize_result.jac.shape) * s[0]
    s = s[s > threshold]
    VT = VT[:s.size]
    pcov = np.dot(VT.T / s**2, VT)
    pcov = pcov * s_sq
           

    if pcov is None:
        # indeterminate covariance
        print('Failed estimating pcov1')
        pcov = np.zeros((len(pcov), len(pcov)), dtype=float)
        pcov.fill(np.inf)
        
    return pcov
    
def fit_rsh(V: np.ndarray, J: np.ndarray):
    all_tol = np.finfo(np.float64).eps
    res = optimize.least_squares(fobj,[1E-8, 1E-3], 
                                 jac=jac,
                                 bounds=([1E-10, 1E-5],[1E50, 1E50]),
                                 args=(J, V),
                                 xtol=all_tol,
                                 ftol=all_tol,
                                 gtol=all_tol,
#                                 x_scale='jac',
#                                 loss='soft_l1', f_scale=0.1,
#                                 loss='cauchy', f_scale=0.1,
                                 max_nfev=len(V)*1000,
                                 verbose=0)
    return res
    

if __name__ == '__main__':
    mpl.rcParams.update(defaultPlotStyle)
    analysis = pia.Analysis(folder_path=base_folder)
    
        
    if platform.system() == 'Windows':
        base_folder = r'\\?\\' + os.path.abspath(base_folder)
    
    results_folder = os.path.join(base_folder, output_folder)
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    
    
    efficiency_data = analysis.batch_efficiency(csv_index=csv_index)
    # How many points
    n_total = len(efficiency_data)
    
    time_s = efficiency_data['time (s)']
    time_h = time_s / 3600
    df = pd.read_csv(filepath_or_buffer=os.path.join(base_folder, 'jv_plots', csv_index))
   
    
    cmap = mpl.cm.get_cmap('viridis_r')
    normalize = mpl.colors.Normalize(vmin=np.amin(time_h), vmax=np.amax(time_h))
    jv_colors = [cmap(normalize(t)) for t in time_h]
    scalar_maps = mpl.cm.ScalarMappable(cmap=cmap, norm=normalize)
    
    with h5py.File(na_file, 'r') as hf:
        group_si = hf['L2']
        x = np.array(group_si['x'])
        x = x - x[0]
       
        

    #plt.ioff()
    fig = plt.figure()
    fig.set_size_inches(6.0, 3.0, forward=True)
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    gs0 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, width_ratios=[1])
    gs00 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, subplot_spec=gs0[0])
    ax1 = fig.add_subplot(gs00[0, 0])
    ax2 = fig.add_subplot(gs00[0, 1])
    rsh_data = np.empty(len(efficiency_data), dtype=rsh_dtype)

    voc_max = np.amax(efficiency_data['voc (V)'])
    with h5py.File(na_file, 'r') as hf:
        group_si = hf['L2']
        profiles = group_si['concentration']
        for i, row in df.iterrows():
            ct_ds = 'ct_{0:d}'.format(row['index'])
            t = float(row['time (s)']) / 3600
            c = np.array(profiles[ct_ds])
            jv = analysis.read_jv(h5_filename=os.path.join(base_folder, 'jv_plots', row['filename']))
            jv = jv[jv['current (mA/cm2)'] >= 0]
            ax2.plot(jv['voltage (V)'], jv['current (mA/cm2)'], '-', color=cmap(normalize(t)))
            ax1.plot(x, c, color=jv_colors[i])
            # Fit Rsh from the interval between 0 and 0.1 V (assume linear behavior)
            idx_fit = jv['voltage (V)'] <= 0.15
            res = fit_rsh(V=jv['voltage (V)'][idx_fit], J=jv['current (mA/cm2)'][idx_fit]/1000)
            popt = res.x
#            print('popt = {0}'.format(popt))
            pcov = find_pcov(res)
            ci = cf.confint(n=len(jv['voltage (V)'][idx_fit]), pars=popt, pcov=pcov, confidence=0.95)
#            ci = np.power(10, ci)
            rsh_data[i] = (time_h[i], popt[0], ci[0,0], ci[0,1])
            

    ax2.plot(efficiency_data['v_mpp (V)'], efficiency_data['j_mpp (mA/cm2)'],
                             'o', color='r', ms=5, alpha=0.5)
    
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.03)
    cbar = fig.colorbar(scalar_maps, cax=cax)
    cbar.set_label('Time (hr)', rotation=90)

    ax1.set_xlabel('Depth (um)')
    ax1.set_ylabel('[Na] (cm$^{-3}$)')
    

    ax2.set_xlabel('Bias (V)')
    ax2.set_ylabel('J (mA/cm$^{2}$)')
    
    ax1.set_yscale('log')
    ax1.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=6))
    ax1.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=60, subs=np.arange(2, 10) * .1))
    ax1.set_xlim(0, np.amax(x))
    ax1.set_ylim(bottom=1E14, top=1E18)
    ax2.set_xlim(0, voc_max)

    ax1.xaxis.set_major_formatter(xfmt)
    ax1.xaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    
    ax2.xaxis.set_major_formatter(xfmt)
    ax2.xaxis.set_major_locator(mticker.MaxNLocator(4, prune=None))
    ax2.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
   
    ax2.yaxis.set_major_formatter(xfmt)
    ax2.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    ax2.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))

    ax1.set_title('Concentration Profile')
    ax2.set_title('Illuminated JV')
    
    

    plt.tight_layout()
    plt.show()
    
    fig.savefig(os.path.join(results_folder, 'illuminated_jv.png'), dpi=600)
    
    
    # Plot the efficiency data
    fig = plt.figure()
    fig.set_size_inches(5.0, 7.0, forward=True)
    fig.subplots_adjust(hspace=0.1, wspace=0.35)
    gs0 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, width_ratios=[1])
    gs00 = gridspec.GridSpecFromSubplotSpec(nrows=4, ncols=1, subplot_spec=gs0[0])
    ax1 = fig.add_subplot(gs00[0, 0])
    ax2 = fig.add_subplot(gs00[1, 0])
    ax3 = fig.add_subplot(gs00[2, 0])
    ax4 = fig.add_subplot(gs00[3, 0])
    
    p0 = efficiency_data['pd_mpp (mW/cm2)'][0]
    nu0 = efficiency_data['efficiency'][0]
    
    ax1.plot(time_h, 100*efficiency_data['pd_mpp (mW/cm2)']/p0, color='C0')
    ax2.fill_between(time_h, rsh_data['Rsh lci'], rsh_data['Rsh uci'], color=lighten_color(mpl.colors.to_rgb('C1')))
    ax2.plot(time_h, rsh_data['Rsh (Ohms cm2)'], color='C1')
    ax3.plot(time_h, 100*efficiency_data['efficiency']/nu0, color='C2')
    ax4.plot(time_h, 1000*efficiency_data['voc (V)'], color='C3')
    
    ax1.set_xlabel('Time (hr)')
    ax4.set_xlabel('Time (hr)')
    
    ax1.set_ylabel('$P_{\\mathrm{mpp}}/P_{\\mathrm{mpp}}(t=0)$')
    ax2.set_ylabel('$R_{\mathrm{sh}}$ ($\Omega$ cm$^{2}$)')
    ax3.set_ylabel('$\\eta/\\eta(t=0)$')
    ax4.set_ylabel('V$_{\\mathregular{oc}}$ (V)')
    
    ax1.set_xlim(0, np.amax(time_h))
    ax2.set_xlim(0, np.amax(time_h))
    ax3.set_xlim(0, np.amax(time_h))
    
    ax1.xaxis.set_label_position('top')
    
    ax1.tick_params(labelbottom=False, top=True, right=True, which='both', labeltop=True)
    ax2.tick_params(labelbottom=False, top=False, right=True, which='both')
    ax3.tick_params(labelbottom=False, top=False, right=True, which='both')
    ax4.tick_params(labelbottom=True, top=False, right=True, which='both')
    
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
   
#    ax2.yaxis.set_major_formatter(xfmt)
#    ax2.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
#    ax2.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    ax3.xaxis.set_major_formatter(xfmt)
    ax3.xaxis.set_major_locator(mticker.MaxNLocator(12, prune=None))
    ax3.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
   
    ax3.yaxis.set_major_formatter(xfmt)
    ax3.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    ax3.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    ax4.xaxis.set_major_formatter(xfmt)
    ax4.xaxis.set_major_locator(mticker.MaxNLocator(12, prune=None))
    ax4.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
   
    ax4.yaxis.set_major_formatter(xfmt)
    ax4.yaxis.set_major_locator(mticker.MaxNLocator(5, prune=None))
    ax4.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    fig.savefig(os.path.join(results_folder, 'efficiency_plot.png'), dpi=600)
    rsh_df = pd.DataFrame(data=rsh_data)
    rsh_df.to_csv(path_or_buf=os.path.join(results_folder, 'rsh_data.csv'), index=False)
    
    
    fig = plt.figure()
    fig.set_size_inches(4.0, 3.7, forward=True)
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    gs0 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, width_ratios=[1])
    gs00 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=1, subplot_spec=gs0[0])
    ax1 = fig.add_subplot(gs00[0, 0])
    for i, row in df.iterrows():
        t = float(row['time (s)']) / 3600
        jv = analysis.read_jv(h5_filename=os.path.join(base_folder, 'jv_plots', row['filename']))
        jv = jv[jv['current (mA/cm2)'] >= 0]
        # Fit Rsh from the interval between 0 and 0.1 V (assume linear behavior)
        idx_fit = jv['voltage (V)'] <= 0.15
        res = fit_rsh(V=jv['voltage (V)'][idx_fit], J=jv['current (mA/cm2)'][idx_fit])
        popt = res.x
#        print('popt = {0}'.format(popt))
        pcov = find_pcov(res)
#        ci = cf.confint(n=len(jv['voltage (V)'][idx_fit]), pars=popt, pcov=pcov, confidence=0.95)
#            ci = np.power(10, ci)
        ax1.plot(jv['current (mA/cm2)'][idx_fit], jv['voltage (V)'][idx_fit], 'o', fillstyle='none', color=cmap(normalize(t)))
        ax1.plot(jv['current (mA/cm2)'][idx_fit], -popt[0]*jv['current (mA/cm2)'][idx_fit] + popt[1], '-', color=cmap(normalize(t)))
        
    ax1.set_xlabel('J (mA/cm$^2$)')
    ax1.set_ylabel('Bias (V)')
    
    fig = plt.figure()
    fig.set_size_inches(5.0, 4.5, forward=True)
    fig.subplots_adjust(hspace=0.1, wspace=0.35)
    gs0 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig, width_ratios=[1])
    gs00 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=1, subplot_spec=gs0[0])
    ax1 = fig.add_subplot(gs00[0, 0])
    ax2 = fig.add_subplot(gs00[1, 0])
    
    
    ax1.plot(time_h, 100*efficiency_data['pd_mpp (mW/cm2)']/p0, color='C0')
    ax2.fill_between(time_h, rsh_data['Rsh lci'], rsh_data['Rsh uci'], color=lighten_color(mpl.colors.to_rgb('C1')))
    ax2.plot(time_h, rsh_data['Rsh (Ohms cm2)'], color='C1')
    
    ax2.set_xlabel('Time (hr)')
    
    ax1.set_ylabel('$P_{\\mathrm{mpp}}/P_{\\mathrm{mpp}}(t=0)$')
    ax2.set_ylabel('$R_{\mathrm{sh}}$ ($\Omega$ cm$^{2}$)')

    
    ax1.set_xlim(0, np.amax(time_h))
    ax2.set_xlim(0, np.amax(time_h))
    
    ax1.xaxis.set_label_position('top')
    
    ax1.tick_params(labelbottom=False, top=True, right=True, which='both')
    ax2.tick_params(labelbottom=True, bottom=True, top=False, right=True, which='both')
    
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
   

    
    fig.savefig(os.path.join(results_folder, 'pmax_plot.png'), dpi=600)
    
        