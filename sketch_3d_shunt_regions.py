import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.animation as manimation
import pid.utils as utils
import matplotlib as mpl
import numpy as np

# Define the output filename
filetag = 'shunt_3d'


def update_line(n, line, ax_, elev=30):
    # ax_: Axes3D = line[0]
    azim = n % 360
    print('Rotating view: elev = {0:.1f}, azim = {1:.1f}'.format(30, azim))
    ax_.view_init(elev=elev, azim=azim)
    return line


if __name__ == "__main__":
    # Define the depth of the shunt in um
    shunt_depth = 1.
    # define the number of regions
    n_regions = 20
    # Get the coordinates of the regions (centered at x=0, y=0, z=0)
    regions = utils.get_triangle_regions(shunt_depth=shunt_depth, n_regions=n_regions)
    # create the figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # create a color map
    cmap = mpl.cm.get_cmap('viridis_r')
    normalize = mpl.colors.Normalize(vmin=0, vmax=n_regions)
    sigma_colors = [cmap(normalize(t)) for t in range(0, n_regions)]
    # create the plot
    line = []
    for i, reg in enumerate(regions):
        vertex = [[p for p in reg['points'][:-1]]]
        center = reg['center']
        polygon = Poly3DCollection(vertex)
        polygon.set_facecolor(sigma_colors[i])
        polygon.set_edgecolor('k')
        ph1 = ax.add_collection3d(polygon)
        ph2 = ax.scatter(center[0], center[1], center[2], color='k', marker='+')
        # line.append(ph1)
        line.append(ph2)

    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_zlim(-1, 0.0)
    ax.set_xlabel('x (um)')
    ax.set_ylabel('y (um)')
    ax.set_zlabel('z (um)')
    ax.set_title('Metallic shunt')
    ax.view_init(elev=30, azim=0)
    plt.tight_layout()

    # line = [ax]
    print(line)
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Sketch of 3D shunt geometry', artist='Matplotlib',
                    comment='author: Erick Martinez Loran')
    writer = FFMpegWriter(fps=30, metadata=metadata, extra_args=['-vcodec', 'libx264'])

    ani = manimation.FuncAnimation(fig=fig, func=update_line, frames=np.arange(0., 181), interval=10, blit=True,
                                   repeat=False, fargs=(line, ax))

    # plt.show()

    ani.save(filetag + '.mp4', writer=writer, dpi=300)

