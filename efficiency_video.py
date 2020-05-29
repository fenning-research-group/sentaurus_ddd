# -*- coding: utf-8 -*-
"""
Created on Sun May  3 17:58:55 2020

@author: Erick
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.animation as manimation
from matplotlib.ticker import ScalarFormatter
import h5py
import os
import pandas as pd
import platform
import matplotlib.gridspec as gridspec
import re
import pid.analysis as pia

#dataPath = r'G:\My Drive\Research\PVRD1\FENICS\SUPG_TRBDF2\simulations\results_two_layers\pnp\SWEEP_D,H'
#dataFile = r'two_layers_D1=4E-16cm2ps_D2=1E-14cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-12_m1.0e+00_pnp.h5'
#results_path = 'videos'

# radius = 50E-4 # cm
area = np.pi*(5)*1E-6


base_folder = r'C:\Users\Erick\PycharmProjects\sentaurus_ddd\results\time_dependence\20200428_sc=1E2-D2=1E-14_shuntDepth=1.0-0'
output_folder = r'analysis_plots'
csv_index = r'file_index.csv'
na_file = r'C:\Users\Erick\PycharmProjects\sentaurus_ddd\results\time_dependence\20200428_sc=1E2-D2=1E-14_shuntDepth=1.0-0\two_layers_D1=4E-16cm2ps_D2=1E-14cm2ps_Cs1E+20cm3_T85_time96hr_h1.0e-10_m1.0e+00_pnp.h5'
xfmt = ScalarFormatter(useMathText=True)
xfmt.set_powerlimits((-3, 3))


max_time = 40 # hr

plot_vfb = False

plotStyle = {'font.size': 12,
                 'font.family': 'Arial',
                 'font.weight': 'regular',
                'legend.fontsize': 12,
                'mathtext.fontset': 'custom',
                'mathtext.rm': 'Times New Roman',
                'mathtext.it': 'Times New Roman:italic',#'Arial:italic',
                'mathtext.cal': 'Times New Roman:italic',#'Arial:italic',
                'mathtext.bf': 'Times New Roman:bold',#'Arial:bold',
                'xtick.direction' : 'in',
                'ytick.direction' : 'in',
                'xtick.major.size' : 4.0,
                'xtick.major.width' : 1.5,
                'ytick.major.size' : 4.0,
                'ytick.major.width' : 1.5,
                'xtick.minor.size' : 2.5,
                'xtick.minor.width' : 0.7,
                'ytick.minor.size' : 2.5,
                'ytick.minor.width' : 0.7,
                'lines.linewidth'   : 3,
                'lines.markersize'  : 10,
                'lines.markeredgewidth'  : 1.0,
                'axes.labelpad'  : 6.0,
                'axes.labelsize' : 12,
                'axes.labelweight' : 'regular',
                'axes.titlesize' : 12,
                'axes.titleweight' : 'bold',
                'axes.titlepad' : 8,
                'figure.titleweight' : 'bold',
                'figure.dpi': 100}

q_red   = 1.6021766208 # x 1E-19 C
e0_red  = 8.854187817620389 # x 1E-12 C^2 / J m


def conductivity_model(concentration: np.ndarray, activated_na: float = 1, 
                       segregation_coefficient: float = 120) -> np.ndarray:
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
    cshunt = concentration * segregation_coefficient * activated_na

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

def latexFormat(x,digits=2):
    fmt_dgts = '%%.%df' % digits
    fmt_in   = '%%.%dE' % digits
    x_str = fmt_in % x
    x_sci = (np.array(x_str.split('E'))).astype(np.float)
    if digits == 0:
        return r'$\mathregular{10^{%d}}$' % x_sci[1]
    else:
        ltx_str = fmt_dgts % x_sci[0]
        ltx_str += r'$\mathregular{\times 10^{%d}}$' % x_sci[1]
        return ltx_str
    
def tau_c(D: float,E: float,L: float,T: float) -> float:
    '''
    tau_c estimates the charactersitic constant for the Nernst-Planck
    equation in the low concentration approximation
    
    tau_c = 2D/(mu^2 E^2) + L/(mu E)*[1 ± 2*(kT/(q E L))^(1/2)]
    
    Since mu = qD/kT
    
    tau_c = (2/D) X^2 + (l/D) X * [1 ± 2*(X/L)^(1/2)],
    
    with X = kT/qE
    
    Parameters:
    -----------
    D: float
        The diffusion coefficient in cm^2/s
    E: float
        The electric field in MV/cm = 1E6 V/cm
    L: float
        The distance in cm
    T: float
        The temperature in °C
    
    Returns:
    --------
    float
        The characteristic time in s
    '''
    TK      = T + 273.15
    kB_red  = 8.6173303 # 1E-5 eV/K
    
    kTq_red = kB_red*TK # x 1E-5 V
    
    x       =  kTq_red/E  # x 1E-5 V x 1E-6 cm/V = 1E-11 cm
    
    tau1    = 1.0E-11*x*(L/D) + (2.0E-22*np.power(x,2.0)/D)*(1.0 - np.sqrt(1.0+1E-11*(L/x)))
    tau2    = 1.0E-11*x*(L/D) + (2.0E-22*np.power(x,2.0)/D)*(1.0 + np.sqrt(1.0+1E-11*(L/x)))
#    
#    E = 1E6*E
#    
#    kTq = constants.value('Boltzmann constant in eV/K')*TK
#    mu  = D/kTq
#    tau1     = 2.0*D/(np.square(mu*E)) \
#                + L/(mu*E)*(1.0-2.0*np.sqrt(kTq/(E*L)))
#    tau2     = 2.0*D/(np.square(mu*E)) \
#                + L/(mu*E)*(1.0+2.0*np.sqrt(kTq/(E*L)))
#    return tau2
    
    if tau1 <= 0:
        return tau2
    elif tau2 <= 0:
        return tau1
    else:
        return min(tau1, tau2)
    


def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \times 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str


def update_line(i, line, h5file, x1, x2, time_s_, tau, sigma, nu, t_idx, nu0_):
    with h5py.File(h5file, 'r') as hf:
        grp_sinx    = hf['/L1']  
        grp_si      = hf['/L2']       
        ct_ds = 'ct_{0:d}'.format(t_idx[i])
        c1 = np.array(grp_sinx['concentration'][ct_ds])
        c2 = np.array(grp_si['concentration'][ct_ds])
        t = time_s_[i]
        
#        print('i={0}, nu = {1}'.format(i,nu[0:i]*100))
        
        line[0].set_data(x1,c1)
        line[1].set_data(x2,c2)
        if i>0:
            line[2].set_data(np.array(time_s_[0:i]/3600), np.array(sigma[0:i]))
            line[3].set_data(np.array(time_s_[0:i]/3600), np.array(nu[0:i]*100/nu0_))
            
            jv = analysis.read_jv(h5_filename=os.path.join(base_folder, 'jv_plots', df.iloc[i]['filename']))
            line[4].set_data(np.array(jv['voltage (V)']), np.array(jv['current (mA/cm2)']))
            
        line[5].set_text(r'{0:.1f} hr'.format(t/3600))
        
    
#    if t >= tau:
#        line[0].set_color('r')
#        line[1].set_color('r')
#    
    

    print('Updating time step {0}/{1}'.format(i, len(time_s_)))
    return line


if __name__ == '__main__':
        
    results_folder = os.path.join(base_folder, output_folder)
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    
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
    
    # Get the indices of the time steps
    t_index = df['index'].astype(int)
    conductivities = np.empty(n_total, dtype=float)
    for i, tix in enumerate(t_index):
        with h5py.File(na_file, 'r') as hf:
            ct_ds = 'ct_{0:d}'.format(tix)
            grp_si      = hf['/L2']   
            c = np.array(grp_si['concentration'][ct_ds])
            c = c[(np.abs(c - c[0]/100)).argmin()]
            conductivities[i] = conductivity_model(concentration=c)
    
    filetag     = os.path.splitext(os.path.basename(na_file))[0]
#    if platform.system() == 'Windows':
#        fullfile = u'\\\?\\' + fullfile
#        output_path = u'\\\?\\' + output_path
        
        
    x1 = None
    x2 = None
    c1 = None
    c2 = None
    
    
    pattern1 = re.compile('h(\d+\.?\d*e[+-]\d{1,2})')
    match1 = pattern1.search(na_file)
    h = float(match1.group(1))
    
    jv0 = analysis.read_jv(h5_filename=os.path.join(base_folder, 'jv_plots', df.iloc[0]['filename']))
        
    with h5py.File(na_file, 'r') as hf:
        
        grp_time    = hf['time']
#        time        = np.array(hf['time'])
    #    vfb = np.zeros_like(time)
        vfb         = np.array(hf['vfb'])
        
        TempC   = grp_time.attrs['temp_c']
        Cs      = grp_time.attrs['Csource']
        Cbulk   = grp_time.attrs['Cbulk']
#        
#        h       = grp_time.attrs['h']
#        m       = grp_time.attrs['m']
        
        
        grp_sinx    = hf['/L1']
        grp_si      = hf['/L2']
        
        er          = grp_sinx.attrs['er']
        E1          = grp_sinx.attrs['electric_field_eff']*er
        D1          = grp_sinx.attrs['D']
        V1          = grp_sinx.attrs['stress_voltage']
        
        D2          = grp_si.attrs['D']
    
           
        
        x1          = np.array(grp_sinx['x'])*1000
        x2          = np.array(grp_si['x'])
        
        L1          = np.amax(x1)
        x1          = x1 - L1
        x2          = x2 - L1/1000
    
        
        tau_c       = tau_c(D1,E1/er,L1*1E-7,TempC) 
        
        
        c1 = np.array(grp_sinx['concentration']['ct_0'])    
        c2 = np.array(grp_si['concentration']['ct_0'])    
    
    
    # Plot style parameters
    mpl.rcParams.update(plotStyle)
    imax = (np.abs(time_s - max_time*3600)).argmin()
    
    fig1 = plt.figure()
    fig1.set_size_inches(6.5,6.5,forward=True)
    fig1.subplots_adjust(hspace=0.35, wspace=0.35)
    gs0 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig1, hspace=0.85)
    gs00 = gridspec.GridSpecFromSubplotSpec(nrows=1, ncols=2, width_ratios=[0.5, 1.0],
                                            subplot_spec = gs0[0], wspace=0)
    gs01 = gridspec.GridSpecFromSubplotSpec(nrows=2, ncols=2, width_ratios=[0.5, 1.0],
                                            subplot_spec = gs0[1::,0], hspace=0.1)
    
    
    
    
    
    gs2 = fig1.add_gridspec(nrows=1, ncols=1)
    ax1 = fig1.add_subplot(gs00[0, 0])
    ax2 = fig1.add_subplot(gs00[0, 1])
    ax3 = fig1.add_subplot(gs01[0, 0])
    ax4 = fig1.add_subplot(gs01[1, 0])
    ax5 = fig1.add_subplot(gs01[:, 1])
    
    
    
    
    
    ax1.set_facecolor((0.89, 0.75, 1.0))
    ax2.set_facecolor((0.82, 0.83, 1.0))
    
    ph1,  = ax1.plot(x1, c1, color='C0')
    ph2,  = ax2.plot(x2, c2, color='C0')
    ph3,  = ax3.plot([], [], color='C1')
    ph4,  = ax4.plot([], [], color='C2')
    ph5,  = ax5.plot(jv0['voltage (V)'], jv0['current (mA/cm2)'], color='C3')
    
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax2.set_yscale('log')
    
    nu0 = efficiency_data['efficiency'][0]
    
    ax1.set_xlim([np.amin(x1),np.amax(x1)])
    ax2.set_xlim([np.amin(x2),np.amax(x2)])
    ax1.set_ylim([1E12,1E22])
    ax2.set_ylim([1E12,1E22])
    
    ax3.set_xlim(0, np.amax(time_h))
    ax3.set_ylim(bottom=np.amin(conductivities)*0.2, top=np.amax(conductivities)*1.5)
    ax4.set_xlim(0, np.amax(time_h))
    ax4.set_ylim(bottom=np.amin(efficiency_data['efficiency'])*100/nu0, top=100)
    
    ax5.set_xlim(np.amin(jv0['voltage (V)']), np.amax(jv0['voltage (V)']))
    ax5.set_ylim(bottom=0.0)
    
    ax1.tick_params(labelbottom=False, bottom=True, top=False, right=False, which='left')
    ax2.tick_params(labelbottom=True, labelleft=False, left=False, bottom=True, top=False, right=True, which='right')
        
    ax3.tick_params(labelbottom=False, top=False, right=True, which='both')
    ax4.tick_params(labelbottom=True, top=False, right=True, which='both')
    
    
    ax1.set_ylabel("[Na$^{+}$] (cm$^{-3}$)")
    ax1.set_xlabel("Depth (nm)")
    
    ax2.set_xlabel("Depth ($\mu$m)")
    
    
    ax3.set_ylabel("$\\sigma$ ([Na]) (S/cm)")
    ax4.set_ylabel("$\\eta (t) $/$\\eta_{t=0}$ (%)")
    
    ax5.set_xlabel('Bias (V)')
    ax5.set_ylabel('J (mA/cm$^2$)')
    ax5.set_title('Illuminated JV')
    
    
    
    ax4.set_xlabel("Time (hr)")
    
    
    time_txt = ax5.text(0.95, 0.95, '0.0', 
                        horizontalalignment='right',
                        verticalalignment='top', 
                        transform=ax5.transAxes)
    
    params1_str = '$C_s$ = %s cm$^{-3}$\nD = %s cm$^2$/s' % (latexFormat(Cs, digits=1),latexFormat(D1, digits=1))
    params1_txt = ax1.text(0.05, 0.95, params1_str, 
                        horizontalalignment='left',
                        verticalalignment='top', 
                        transform=ax1.transAxes,
                        fontsize=11,
                        color='b')
    
    h_exp  = int(("%e" % h).split('e')[1])
    params2_str = 'D = %s cm$^2$/s\n$h$ = 10$^{%d}$ cm/2' % (latexFormat(D2, digits=1), h_exp)
    params2_txt = ax2.text(0.05, 0.95, params2_str, 
                        horizontalalignment='left',
                        verticalalignment='top', 
                        transform=ax2.transAxes,
                        fontsize=11,
                        color='b')
    
    
    
    ax1.set_title('SiN$_{\mathregular{x}}$, %s MV/cm' % (latex_float(E1)))
    ax2.set_title('Si, E = 0')
    
    xfmt = ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((-4,4))
    locmaj1 = mpl.ticker.LogLocator(base=10.0, numticks=5) 
    locmin1 = mpl.ticker.LogLocator(base=10.0, numticks=50, subs=np.arange(2, 10) * .1) 
    
    
        
    ax1.xaxis.set_major_locator(mticker.MaxNLocator(5,prune='upper'))
    ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    
       
    ax1.yaxis.set_major_locator(locmaj1)
    ax1.yaxis.set_minor_locator(locmin1)
    
    # Axis 2
    locmaj2 = mpl.ticker.LogLocator(base=10.0, numticks=5)
    locmin2 = mpl.ticker.LogLocator(base=10.0, numticks=50, subs=np.arange(2, 10) * .1)
    ax2.yaxis.set_ticks_position('right')
    
    ax2.xaxis.set_major_locator(mticker.MaxNLocator(5,prune=None))
    ax2.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    
       
    ax2.yaxis.set_major_locator(locmaj2)
    ax2.yaxis.set_minor_locator(locmin2)
#    ax1.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    

    ax3.xaxis.set_major_locator(mticker.MaxNLocator(5,prune=None))
    ax3.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    ax3.yaxis.set_major_locator(mticker.MaxNLocator(6,prune='lower'))
    ax3.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    
    ax4.xaxis.set_major_locator(mticker.MaxNLocator(5,prune=None))
    ax4.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    ax4.yaxis.set_major_locator(mticker.MaxNLocator(4,prune='lower'))
    ax4.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    ax5.xaxis.set_major_locator(mticker.MaxNLocator(4,prune=None))
    ax5.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    ax5.yaxis.set_major_locator(mticker.MaxNLocator(7,prune='lower'))
    ax5.yaxis.set_minor_locator(mticker.AutoMinorLocator(2))

    
#    plt.tight_layout()
#    plt.show()
    
    line = [ph1, ph2, ph3, ph4, ph5, time_txt]
    
    
    
    
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Simulated SIMS profile', artist='Matplotlib',
                    comment='Time dependent profile')
    writer = FFMpegWriter(fps=1,metadata=metadata, extra_args=['-vcodec', 'libx264'])
    ani = manimation.FuncAnimation(fig1, update_line, blit=True, interval=200,
        repeat=False,frames=np.arange(0,n_total),
        fargs=(line, na_file, x1, x2, np.array(time_s), tau_c, conductivities, 
               np.array(efficiency_data['efficiency']), np.array(t_index), nu0))
    
    
    ft = os.path.join(results_folder, 'device_simulation.mp4')
#    if platform.system() == 'Windows':
#        ft = u'\\\?\\' + ft
    ani.save(ft, writer=writer,dpi=300)
    
    