# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:24:24 2019

@author: Erick
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.animation as manimation
from matplotlib.ticker import ScalarFormatter
import numdiffusion as nmd
import h5py
import os

dataPath = r'G:\My Drive\Research\PVRD1\FENICS\SUPG_TRBDF2\PNP\twolayers\simulations'
dataFile = r'two_layer_stack_h_1E-4m_10_results.h5'


plotStyle = {'font.size': 18,
                 'font.family': 'Arial',
                 'font.weight': 'regular',
                'legend.fontsize': 16,
                'mathtext.fontset': 'custom',
                'mathtext.rm': 'Arial',
                'mathtext.it': 'Arial:italic',#'Arial:italic',
                'mathtext.cal': 'Arial:italic',#'Arial:italic',
                'mathtext.bf': 'Arial:bold',#'Arial:bold',
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
                'axes.labelsize' : 18,
                'axes.labelweight' : 'regular',
                'axes.titlesize' : 20,
                'axes.titleweight' : 'bold',
                'axes.titlepad' : 8,
                'figure.titleweight' : 'bold',
                'figure.dpi': 300}


def update_line(i,line,hf,x1,x2,time_s):
    grp_sinx    = hf['sinx']
    grp_si      = hf['si']
    
    L1          = np.amax(x1)
    
    ct_ds = 'ct_{0:d}'.format(i)
    c1 = np.array(grp_sinx['concentration'][ct_ds])
    c2 = np.array(grp_si['concentration'][ct_ds])
    t = time_s[i]
    line[0].set_data(x1,c1)
    line[1].set_data(x2-L1,c2)
    line[2].set_text('{0}'.format(nmd.formatTimeStr(t)))
    print('Updating time step {0}/{1}.'.format(i,len(time_s)))
    return line

if __name__ == '__main__':
    filetag = ''.join((dataFile.split('.'))[0:-1])
    hf          = h5py.File(dataFile, 'r')
    hfKeys      = list(hf.keys())
    
    grp_time    = hf['time']
    time        =  np.array(hf['time'])
    
    TempC   = grp_time.attrs['temp_c']
    Cs      = grp_time.attrs['Csource']
    Cbulk   = grp_time.attrs['Cbulk']
    h       = grp_time.attrs['h']
    m       = grp_time.attrs['m']
    
    grp_sinx    = hf['sinx']
    grp_si      = hf['si']
    
    E1          = grp_sinx['x'].attrs['electric_field']
    D1          = grp_sinx['x'].attrs['D']
    E2          = grp_si['x'].attrs['electric_field']
    D2          = grp_si['x'].attrs['D']
    
    x1          = np.array(hf['sinx']['x'])*1000
    x2          = np.array(hf['si']['x'])*1000
    
    L1          = np.amax(x1)
    L2          = np.amax(x2)
    
    c1 = np.array(grp_sinx['concentration']['ct_0'])
    c2 = np.array(grp_si['concentration']['ct_0'])
    
    # Plot style parameters
    mpl.rcParams.update(nmd.defaultPlotStyle)
    fig1 = plt.figure(figsize=(800,400))
    fig1.set_size_inches(8.0,4.0,forward=True)
    ax1 = plt.subplot2grid((1,2), (0,0),fig=fig1)
    ax2 = plt.subplot2grid((1,2), (0,1),fig=fig1)
    
    
    
        
     # Plot solution
    #ph1,  = ax1.semilogy(x_1,y_1,marker='o',color='C0',fillstyle='none')
    #ph2,  = ax2.semilogy(x_2,y_2,marker='o',color='C1',fillstyle='none')
    #ph3,  = ax3.semilogy(x_3,y_3,marker='o',color='C2',fillstyle='none')
    
    ph1,  = ax1.semilogy(x1,c1,color='C1')
    ph2,  = ax2.semilogy((x2-L1),c2,color='C2')
    
    ax1.set_xlim([0,L1])
    ax1.set_ylim([1E12,1E22])
    
    ax2.set_xlim([0,400])
    ax2.set_ylim([1E12,1E22])
    
    
    ax1.set_ylabel("[Na$^{+}$] (cm$^{-3}$)")
    ax1.set_xlabel("Depth (nm)")
    ax2.set_xlabel("Depth (nm)")
    
    time_txt = ax1.text(0.05, 0.05, nmd.formatTimeStr(0.0), 
                        horizontalalignment='left',
                        verticalalignment='bottom', 
                        transform=ax1.transAxes)
    
    h0_exp  = int(("%e" % h).split('e')[1])
    
    h0_str = "h$_{\mathregular{0}}$ = 10$^{%d}$ um/s, m$_{\mathregular{0}}$ = 1\nD = %s cm$^2$/s" % (h0_exp,nmd.latexFormat(D1,2))
    h0_txt = ax1.text(0.05, 0.95, h0_str, 
                        horizontalalignment='left',
                        verticalalignment='top', 
                        transform=ax1.transAxes,
                        fontsize=13,
                        color='C1')
    
    h_exp  = int(("%e" % h).split('e')[1])
    hm_str = "h$_{\mathregular{1}}$ = 10$^{%d}$ um/s, m$_{\mathregular{1}}$ = %g\nD = %s cm$^2$/s" % (h_exp,m,nmd.latexFormat(D2,2))
    hm_txt = ax2.text(0.05, 0.95, hm_str, 
                        horizontalalignment='left',
                        verticalalignment='top', 
                        transform=ax2.transAxes,
                        fontsize=13,
                        color='C2')
    
    ax1.set_title('SiN$_{\mathregular{x}}$, E = %s MV/cm' % (nmd.latex_float(E1)))
    ax2.set_title(r'Si, E = 0')
    
    xfmt = ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((-4,4))
    locmaj = mpl.ticker.LogLocator(base=10.0,numticks=8) 
    locmin = mpl.ticker.LogLocator(base=10.0,numticks=11, subs=np.arange(0.0,1.1,0.1)) 
    
    
    
    ax1.yaxis.set_ticks_position('both')
    ax2.yaxis.set_ticks_position('both')
    
    ax1.xaxis.set_major_locator(mticker.MaxNLocator(7,prune=None))
    ax1.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    
    ax2.xaxis.set_major_locator(mticker.MaxNLocator(5,prune=None))
    ax2.xaxis.set_minor_locator(mticker.AutoMinorLocator(2))
    ax2.xaxis.set_major_formatter(xfmt)
    
    
       
    ax1.yaxis.set_major_locator(locmaj)
    ax1.yaxis.set_minor_locator(locmin)
#    ax1.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    
    ax2.yaxis.set_major_locator(locmaj)
    ax2.yaxis.set_minor_locator(locmin)
#    ax2.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    
    
    
    
    
    
    
    
#    ax1.set_title(r'Concentration profile')
#    ax2.set_title(r'Potential profile')
    
    
    
    
#    ax2.xaxis.set_major_locator(mticker.MaxNLocator(5,prune=None))
#    ax2.xaxis.set_minor_locator(mticker.AutoMinorLocator(5))
#    ax2.yaxis.set_major_formatter(xfmt)
#    ax2.yaxis.set_ticks_position('both')
    
    plt.tight_layout()
    plt.show()
    
    line = [ph1,ph2,time_txt]
    
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Simulated SIMS profile', artist='Matplotlib',
                    comment='Time dependent profile')
    writer = FFMpegWriter(fps=10,metadata=metadata, extra_args=['-vcodec', 'libx264'])
    ani = manimation.FuncAnimation(fig1, update_line, blit=True, interval=200,
        repeat=False,frames=np.arange(0,len(time)),
        fargs=(line,hf,x1,x2,time))
    
    ft = os.path.join(dataPath,filetag)
    ani.save(filetag + '.mp4', writer=writer)
     