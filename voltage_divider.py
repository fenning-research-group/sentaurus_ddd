import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
import platform
import os
import matplotlib.gridspec as gridspec

# Resistivities in Ohm cm
rho_glass_range = [3E12, 1E15]  # [Nagel(a) (soda lime glass), 3M]
# (a) Nagel et al, Quantitative assessment of the local leakage current in PV modules for degradation prediction,
# EUPVSEC 2015
# (b) Naumann et al, On the discrepancy between leakage currents and potential-induced degradation of crystalline
# silicon modules, EUPVSEC 2013
rho_eva_range = [5E13, 1E15]  # [von Gastrow conductivity measurements, GoodFellow]
# See J. Appl. Phys. 49(5), May 1978 (Electrical properties of Si-N films deposited on silicon from
# reactive plasma)
rho_sinx_range = [4E4, 5E19]

# Thicknesses in cm
t_glass = 0.32
t_eva = 45E-4
t_sinx = 75E-7

bias = 1000  # V

output_folder = r'G:\My Drive\Research\PVRD1\Papers\TCAD_draft'


def e_field(layer_index: int, rho: np.ndarray, thickness: np.ndarray, voltage: float = 1000):
    """
    Estimates the electric field in the corresponding layer provided the resistivities and thickness of each of the
    layers and the voltage drop across all the layers

    Parameters
    ----------
    layer_index: int
        The index of the layer on which we estimate the electric field
    rho: np.ndarray
        The resistivities of each of the layers (Ohm cm)
    thickness: np.ndarray
        The thickness of each of the layers (cm)
    voltage: float
        The bias applied to the layer stack

    Returns
    -------
    float
        The electric field in the layer (V / cm)
    """
    return rho[int(layer_index)] * voltage / np.sum(rho * thickness)


def latex_format(number: float, digits: int = 0):
    """
    Formats a number as a latex string in scientific notation

    Parameters
    ----------
    number: float
        The number to format
    digits: int
        The number of digits to round the number to (default = 0)

    Returns
    -------

    """
    number_str = '{0:.9E}'.format(number)
    number_parts = number_str.split('E')
    prefactor = float(number_parts[0])
    exponent = int(number_parts[1])
    format_string = r'{0:.%df} \times 10^{{{1:d}}}' % digits
    return format_string.format(prefactor, exponent)

defaultPlotStyle = {'font.size': 12,
                    'font.family': 'Arial',
                    'font.weight': 'regular',
                    'legend.fontsize': 11,
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


if __name__ == "__main__":
    if platform.system() == 'Windows':
        output_folder = r'\\?\\' + output_folder

    mpl.rcParams.update(defaultPlotStyle)

    thicknesses = np.array([t_glass, t_eva, t_sinx])

    rho_glass = np.logspace(np.log10(rho_glass_range[0]), np.log10(rho_glass_range[1]), 100)
    rho_eva = np.logspace(np.log10(rho_eva_range[0]), np.log10(rho_eva_range[1]), 3)
    rho_sinx = np.logspace(np.log10(rho_sinx_range[0]), np.log10(rho_sinx_range[1]), 100)

    n_glass, n_eva, n_sinx = len(rho_glass), len(rho_eva), len(rho_sinx)

    # Plot the electric field in EVA
    xx, yy = np.meshgrid(rho_glass, rho_sinx)
    e_eva_1 = np.empty((n_glass, n_sinx), dtype=np.float64)
    e_eva_2 = np.empty((n_glass, n_sinx), dtype=np.float64)
    e_eva_3 = np.empty((n_glass, n_sinx), dtype=np.float64)

    fig = plt.figure()
    fig.set_size_inches(7.0, 5.5, forward=True)
    fig.subplots_adjust(hspace=0.5, wspace=0.15)
    gs0 = gridspec.GridSpec(ncols=1, nrows=2, figure=fig)
    gs00 = gridspec.GridSpecFromSubplotSpec(
        nrows=1, ncols=4, subplot_spec=gs0[0], hspace=0.1, width_ratios=[10, 10, 10, 1]
    )
    gs10 = gridspec.GridSpecFromSubplotSpec(
        nrows=1, ncols=4, subplot_spec=gs0[1], hspace=0.1, width_ratios=[10, 10, 10, 1]
    )

    ax00 = fig.add_subplot(gs00[0, 0])
    ax01 = fig.add_subplot(gs00[0, 1])
    ax02 = fig.add_subplot(gs00[0, 2])
    ax03 = fig.add_subplot(gs00[0, 3])

    ax10 = fig.add_subplot(gs10[0, 0])
    ax11 = fig.add_subplot(gs10[0, 1])
    ax12 = fig.add_subplot(gs10[0, 2])
    ax13 = fig.add_subplot(gs10[0, 3])

    ax00.set_xscale('log')
    ax00.set_yscale('log')
    ax01.set_xscale('log')
    ax01.set_yscale('log')
    ax02.set_xscale('log')
    ax02.set_yscale('log')

    ax10.set_xscale('log')
    ax10.set_yscale('log')
    ax11.set_xscale('log')
    ax11.set_yscale('log')
    ax12.set_xscale('log')
    ax12.set_yscale('log')

    ax00.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=3))
    ax00.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=30, subs=np.arange(2, 10) * .1))

    ax01.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=3))
    ax01.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=30, subs=np.arange(2, 10) * .1))

    ax02.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=3))
    ax02.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=30, subs=np.arange(2, 10) * .1))

    ax00.tick_params(
        labelbottom=True, top=True, right=True, which='both', labeltop=False, bottom=True
    )
    ax01.tick_params(
        labelbottom=True, top=True, right=True, which='both', labeltop=False, bottom=True,
        labelleft=False, left=False
    )
    ax02.tick_params(
        labelbottom=True, top=True, right=True, which='both', labeltop=False, bottom=True,
        labelleft=False, left=False
    )

    ##
    ax10.tick_params(
        labelbottom=True, top=True, right=True, which='both', labeltop=False, bottom=True
    )
    ax11.tick_params(
        labelbottom=True, top=True, right=True, which='both', labeltop=False, bottom=True,
        labelleft=False, left=False
    )
    ax12.tick_params(
        labelbottom=True, top=True, right=True, which='both', labeltop=False, bottom=True,
        labelleft=False, left=False
    )

    # First value of rho_eva
    for i, x in enumerate(xx):
        for j, y in enumerate(yy):
            e_eva_1[i, j] = e_field(
                layer_index=1, rho=np.array([rho_glass[i], rho_eva[0], rho_sinx[j]]), thickness=thicknesses,
                voltage=bias
            )

    e_eva_min, e_eva_max = np.amin(e_eva_1), np.amax(e_eva_1)

    # Second value of rho_eva
    for i, x in enumerate(xx):
        for j, y in enumerate(yy):
            e_eva_2[i, j] = e_field(
                layer_index=1, rho=np.array([rho_glass[i], rho_eva[1], rho_sinx[j]]), thickness=thicknesses,
                voltage=bias
            )

    e_eva_min, e_eva_max = min(e_eva_min, np.amin(e_eva_2)), max(e_eva_max, np.amax(e_eva_2))

    # Third value of rho_eva
    for i, x in enumerate(xx):
        for j, y in enumerate(yy):
            e_eva_3[i, j] = e_field(
                layer_index=1, rho=np.array([rho_glass[i], rho_eva[2], rho_sinx[j]]), thickness=thicknesses,
                voltage=bias
            )
    e_eva_min, e_eva_max = min(e_eva_min, np.amin(e_eva_3)), max(e_eva_max, np.amax(e_eva_3))

    normalize1 = mpl.colors.LogNorm(vmin=e_eva_min, vmax=e_eva_max)

    pcm01 = ax00.pcolormesh(xx, yy, e_eva_1, norm=normalize1)
    pcm02 = ax01.pcolormesh(xx, yy, e_eva_2, norm=normalize1)
    pcm03 = ax02.pcolormesh(xx, yy, e_eva_3, norm=normalize1)

    ax00.set_xlabel(r'$\rho_{\mathrm{glass}}$ ($\Omega \cdot \mathrm{cm}$)')
    ax00.set_ylabel(r'$\rho_{\mathrm{SiNx}}$ ($\Omega \cdot \mathrm{cm}$)')
    sci_not = str('{0:.9E}'.format(rho_eva[0])).split('E')
    pf = float(sci_not[0])
    ex = int(sci_not[1])
    ax00.set_title(r'$\rho_{{\mathrm{{EVA}}}} = {0:.0f} \times 10^{{{1:d}}}$ ($\Omega\cdot\mathrm{{cm}}$)'.format(
        pf, ex
    ))

    ax01.set_xlabel(r'$\rho_{\mathrm{glass}}$ ($\Omega \cdot \mathrm{cm}$)')
    sci_not = str('{0:.9E}'.format(rho_eva[1])).split('E')
    pf = float(sci_not[0])
    ex = int(sci_not[1])
    ax01.set_title(r'$\rho_{{\mathrm{{EVA}}}} = {0:.0f} \times 10^{{{1:d}}}$ ($\Omega\cdot\mathrm{{cm}}$)'.format(
        pf, ex
    ))

    ax02.set_xlabel(r'$\rho_{\mathrm{glass}}$ ($\Omega \cdot \mathrm{cm}$)')
    sci_not = str('{0:.9E}'.format(rho_eva[2])).split('E')
    pf = float(sci_not[0])
    ex = int(sci_not[1])
    ax02.set_title(r'$\rho_{{\mathrm{{EVA}}}} = {0:.0f} \times 10^{{{1:d}}}$ ($\Omega\cdot\mathrm{{cm}}$)'.format(
        pf, ex
    ))
    
    ax00.axhline(y=1E11, ls='--', lw=1.0, color='k')
    ax01.axhline(y=1E11, ls='--', lw=1.0, color='k')
    ax02.axhline(y=1E11, ls='--', lw=1.0, color='k')

    ax00.text(
        ax00.get_xlim()[0], 0.7E11, ' Naumann et al.', horizontalalignment='left', verticalalignment='top',
        color='k', fontsize=10
    )

    e_eva_lower = '${0}$ V/cm'.format(latex_format(e_field(
        layer_index=1, rho=np.array([rho_glass[0], rho_eva_range[0], 1E11]), thickness=thicknesses,
        voltage=bias
    ), digits=1))

    e_eva_upper = '${0}$ V/cm'.format(latex_format(e_field(
        layer_index=1, rho=np.array([rho_glass[0], rho_eva_range[1], 1E11]), thickness=thicknesses,
        voltage=bias
    ), digits=1))

    # bbox_props = dict(boxstyle="square,pad=0.15", fc="w", ec="k", lw=0.75)

    ax00.annotate(
        e_eva_lower, xy=(rho_glass[0], 1E11), xycoords='data',
        xytext=(0.2, 0.6), textcoords='axes fraction',
        color='w', fontsize=10,
        # bbox=bbox_props,
        arrowprops=dict(arrowstyle="->",
                        connectionstyle="angle3,angleA=0,angleB=90",
                        color='r'),
    )

    ax02.annotate(
        e_eva_upper, xy=(rho_glass[0], 1E11), xycoords='data',
        xytext=(0.2, 0.6), textcoords='axes fraction',
        color='k', fontsize=10,
        # bbox=bbox_props,
        arrowprops=dict(arrowstyle="->",
                        connectionstyle="angle3,angleA=0,angleB=90",
                        color='r'),
    )

    cbar1 = plt.colorbar(pcm03, pad=0.225, cax=ax03)
    cbar1.set_label(r'$E_{\mathrm{EVA}}$ (V/cm)', rotation=90)

    ax03.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=3))
    ax03.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=30, subs=np.arange(2, 10) * .1))

    del e_eva_1, e_eva_2, e_eva_3

    rho_eva = np.logspace(np.log10(rho_glass_range[0]), np.log10(rho_glass_range[1]), 100)
    rho_glass = np.logspace(np.log10(rho_glass_range[0]), np.log10(rho_glass_range[1]), 3)
    n_eva = len(rho_eva)

    # Plot the electric field in EVA
    xx, yy = np.meshgrid(rho_eva, rho_sinx)
    e_glass_1 = np.empty((n_eva, n_sinx), dtype=np.float64)
    e_glass_2 = np.empty((n_eva, n_sinx), dtype=np.float64)
    e_glass_3 = np.empty((n_eva, n_sinx), dtype=np.float64)

    # First value of rho_glass
    for i, x in enumerate(xx):
        for j, y in enumerate(yy):
            e_glass_1[i, j] = e_field(
                layer_index=1, rho=np.array([rho_glass[0], rho_eva[i], rho_sinx[j]]), thickness=thicknesses,
                voltage=bias
            )

    e_glass_min, e_glass_max = np.amin(e_glass_1), np.amax(e_glass_1)

    # Second value of rho_glass
    for i, x in enumerate(xx):
        for j, y in enumerate(yy):
            e_glass_2[i, j] = e_field(
                layer_index=1, rho=np.array([rho_glass[1], rho_eva[i], rho_sinx[j]]), thickness=thicknesses,
                voltage=bias
            )

    e_glass_min, e_glass_max = min(e_glass_min, np.amin(e_glass_2)), max(e_glass_max, np.amax(e_glass_2))

    # Third value of rho_glass
    for i, x in enumerate(xx):
        for j, y in enumerate(yy):
            e_glass_3[i, j] = e_field(
                layer_index=1, rho=np.array([rho_glass[2], rho_eva[i], rho_sinx[j]]), thickness=thicknesses,
                voltage=bias
            )
    e_glass_min, e_glass_max = min(e_glass_min, np.amin(e_glass_3)), max(e_glass_max, np.amax(e_glass_3))

    normalize2 = mpl.colors.LogNorm(vmin=e_glass_min, vmax=e_glass_max)

    pcm11 = ax10.pcolormesh(xx, yy, e_glass_1, norm=normalize2)
    pcm12 = ax11.pcolormesh(xx, yy, e_glass_2, norm=normalize2)
    pcm13 = ax12.pcolormesh(xx, yy, e_glass_3, norm=normalize2)

    ax10.set_xlabel(r'$\rho_{\mathrm{EVA}}$ ($\Omega \cdot \mathrm{cm}$)')
    ax10.set_ylabel(r'$\rho_{\mathrm{SiNx}}$ ($\Omega \cdot \mathrm{cm}$)')
    sci_not = str('{0:.9E}'.format(rho_glass[0])).split('E')
    pf = float(sci_not[0])
    ex = int(sci_not[1])
    ax10.set_title(r'$\rho_{{\mathrm{{glass}}}} = {0:.0f} \times 10^{{{1:d}}}$ ($\Omega\cdot\mathrm{{cm}}$)'.format(
        pf, ex
    ))

    ax11.set_xlabel(r'$\rho_{\mathrm{EVA}}$ ($\Omega \cdot \mathrm{cm}$)')
    sci_not = str('{0:.9E}'.format(rho_glass[1])).split('E')
    pf = float(sci_not[0])
    ex = int(sci_not[1])
    ax11.set_title(r'$\rho_{{\mathrm{{glass}}}} = {0:.0f} \times 10^{{{1:d}}}$ ($\Omega\cdot\mathrm{{cm}}$)'.format(
        pf, ex
    ))

    ax12.set_xlabel(r'$\rho_{\mathrm{EVA}}$ ($\Omega \cdot \mathrm{cm}$)')
    sci_not = str('{0:.9E}'.format(rho_glass[2])).split('E')
    pf = float(sci_not[0])
    ex = int(sci_not[1])
    ax12.set_title(r'$\rho_{{\mathrm{{glass}}}} = {0:.0f} \times 10^{{{1:d}}}$ ($\Omega\cdot\mathrm{{cm}}$)'.format(
        pf, ex
    ))

    ax10.axhline(y=1E11, ls='--', lw=1.0, color='k')
    ax11.axhline(y=1E11, ls='--', lw=1.0, color='k')
    ax12.axhline(y=1E11, ls='--', lw=1.0, color='k')

    ax10.text(
        ax00.get_xlim()[0], 0.7E11, ' Naumann et al.', horizontalalignment='left', verticalalignment='top',
        color='k', fontsize=10
    )

    e_sinx_upper = '${0}$ V/cm'.format(latex_format(e_field(
        layer_index=2, rho=np.array([rho_glass_range[0], rho_eva[0], 1E11]), thickness=thicknesses,
        voltage=bias
    ), digits=1))

    e_sinx_lower = '${0}$ V/cm'.format(latex_format(e_field(
        layer_index=2, rho=np.array([rho_glass_range[1], rho_eva[0], 1E11]), thickness=thicknesses,
        voltage=bias
    ), digits=1))

    # bbox_props = dict(boxstyle="square,pad=0.15", fc="w", ec="k", lw=0.75)

    ax10.annotate(
        e_sinx_upper, xy=(rho_eva[0], 1E11), xycoords='data',
        xytext=(0.2, 0.6), textcoords='axes fraction',
        color='k', fontsize=10,
        # bbox=bbox_props,
        arrowprops=dict(arrowstyle="->",
                        connectionstyle="angle3,angleA=0,angleB=90",
                        color='r'),
    )

    ax12.annotate(
        e_sinx_lower, xy=(rho_eva[0], 1E11), xycoords='data',
        xytext=(0.2, 0.6), textcoords='axes fraction',
        color='w', fontsize=10,
        # bbox=bbox_props,
        arrowprops=dict(arrowstyle="->",
                        connectionstyle="angle3,angleA=0,angleB=90",
                        color='r'),
    )

    cbar2 = plt.colorbar(pcm13, pad=0.225, cax=ax13)
    cbar2.set_label(r'$E_{\mathrm{SiN_x}}$ (V/cm)', rotation=90)

    ax13.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, numticks=3))
    ax13.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, numticks=30, subs=np.arange(2, 10) * .1))

    del e_glass_1, e_glass_2, e_glass_3

    filetag = 'voltage_divider'

    fig.savefig(os.path.join(output_folder, filetag + '.png'), dpi=600)
    fig.savefig(os.path.join(output_folder, filetag + '.eps'), dpi=600, format='eps')

    plt.tight_layout()
    plt.show()


