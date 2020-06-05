import os
from typing import List
import numpy as np
# import pandas as pd


def files_with_extension(path: str, extension: str) -> List[str]:
    if not os.path.exists(path):
        raise FileNotFoundError('Could not find path: \'{0}\'.'.format(path))
    return [f for f in os.listdir(path) if f.endswith(extension)]


def next_file_name(path, filename) -> str:
    if not os.path.exists(path):
        raise FileNotFoundError('Could not find path: \'{0}\'.'.format(path))
    files = [f for f in os.listdir(path) if f.startswith(filename)]
    n = len(files)
    return '{0}-{1:d}'.format(filename, n)


def make_new_folder(path: str) -> str:
    parent_path = os.path.dirname(os.path.abspath(path))
    if not os.path.exists(parent_path):
        raise FileNotFoundError('Could not find path: \'{0}\'.'.format(parent_path))
    basename = os.path.basename(path)
    folder_name = next_file_name(parent_path, basename)
    full_path = os.path.join(parent_path, folder_name)
    os.makedirs(full_path)
    return full_path


def diffusion_length_spaced(t_max: float, diffusion_coefficient: float, steps: int, t0: float = 0) -> np.ndarray:
    """
    Makes an array of time points spaced according to the time of a uniformly spaced diffusion profile

    t = L^2/D

    Parameters
    ----------
    t_max: float
        The maximum time to simulate (s)
    diffusion_coefficient: float
        The diffusion coefficient in cm^2/s
    steps: int
        The number of steps in the array
    t0: float
        The first time to use (s). Defaults to 0.0


    Returns
    -------
    np.ndarray
        The time points
    """
    Dums = 1E8 * diffusion_coefficient
    max_depth = np.sqrt(Dums * t_max)
    min_depth = np.sqrt(Dums * t0)

    distance = np.linspace(min_depth, max_depth, steps + 1)
    time = np.empty(steps + 1)
    for i, x in enumerate(distance):
        time[i] = (x ** 2) / Dums
    return time


def geometric_series_spaced(max_val: float, min_delta: float, steps: int, reverse: bool = False,
                            **kwargs) -> np.ndarray:
    """
    Produces an array of values spaced according to a geometric series

    S_n = a + a*r + a*r^2 + ... + a*r^(n-2) + a*r^(n-1)

    For which S_n = a*(1-r^n)/(1-r)

    Here, a is the minimum increment (min_delta) and n is the number of steps and r is determined using Newton's method
    
    Parameters
    ----------
    max_val: float
        The maximum value the series will take
    min_delta: float
        The minimum delta value it will take
    steps: int
        The number of values in the array
    reverse: bool
        If true solve for r = 1 / p
    **kwargs: keyword arguments
        n_iterations: int
            The number of Newton iterations to estimate r

    Returns
    -------
    np.ndarray
        A vector with geometrically spaced values
    """
    niter = kwargs.get('n_iterations', 1000)
    debug = kwargs.get('debug', False)
    a = min_delta
    sn = max_val
    n = steps

    if n == 1:
        return np.array([0, max_val])

    # Use Newton's method to determine r

    def f(r_: float) -> float:
        return a * r_ ** n - sn * r_ + sn - a

    def fp(r_: float) -> float:
        return a * n * r_ ** (n - 1) - sn

    r = 2.0
    # Estimate r

    for i in range(niter):
        if fp(r) == 0:  # Avoid zeros
            r += 0.1
        r -= f(r) / fp(r)

    result = np.zeros(n + 1)
    s = 0
    p = 1
    for i in range(n):
        s += a * p
        result[i + 1] = s
        p *= r

    if debug:
        print('r = {0:.3g}'.format(r))

    if reverse:
        new_result = np.zeros_like(result)
        dx = np.diff(result)
        m = len(dx)
        s = 0
        for i in range(m):
            s += dx[m - 1 - i]
            new_result[i + 1] = s
        result = new_result

    return result


def get_indices_at_values(x: np.array, requested_values: np.array) -> np.ndarray:
    """
    Constructs an array of valid indices in the x array corresponding to the requested values

    Parameters
    ----------
    x: np.array
        The array from which the indices will be drawn
    requested_values: np.array

    Returns
    -------
    np.array
        An array with the indices corresponding to the requested values
    """
    result = np.empty(len(requested_values), dtype=int)
    for i, v in enumerate(requested_values):
        result[i] = int((np.abs(v - x)).argmin())
    return result


def get_rectangle_regions(shunt_depth: float = 1, shunt_length: float = 0, n_regions: int = 100, x_center: float = 0,
                          y_center: float = 0):
    """
    Defines the shunt regions for a 3D device

    Parameters
    ----------
    shunt_depth: float
        The depth of the shunt in um
    shunt_length: float
        The length of the rectangular shunt from top (um).
    n_regions: int
        The number of conductivity regions
    x_center: float
        The x position of the base side of the shunt (um)
    y_center: float
        The x position of the base side of the shunt (um)

    Returns
    -------

    """
    # The side of the triangle
    shunt_length = shunt_length if shunt_length > 0 else shunt_depth
    L = np.sqrt(2.) * shunt_depth
    dy = np.sqrt(3.) * L / n_regions / 2.
    z_offset = 0

    data = []
    x_l, x_r = - 0.5 * shunt_length, 0.5 * shunt_length
    isqrt3 = 1.0 / np.sqrt(3.)
    sqrt_23 = np.sqrt(2. / 3.)
    y = 0.0
    point_dtype = np.dtype([('x', 'd'), ('y', 'd'), ('z', 'd')])
    for i in range(n_regions):
        y_next = y + dy

        points = np.empty(5, dtype=point_dtype)
        points[0] = (x_l, y, z_offset)
        points[1] = (x_r, y, z_offset)
        points[2] = (x_r, y_next, z_offset)
        points[3] = (x_l, y_next, z_offset)
        points[4] = (x_l, y, z_offset)
        data.append({
            'region': i,
            'points': points,
            'center': np.array([0.5 * (x_r + x_l), 0.5 * (y_next + y), z_offset])
        })

        y = y_next

    # Rotate the coordinates around the x-axis 54 degrees
    #
    # R =  [ cos(-54.73) -sin(-54.73) ] = [ 1/sqrt(3)    sqrt(2/3) ]
    #      [ sin(-54.73)  sin(-54.73) ]   [ -sqrt(2/3)   1/sqrt(3) ]
    rotation = np.array([
        [1., 0., 0],
        [0., isqrt3, sqrt_23],
        [0., -sqrt_23, isqrt3]
    ])
    for i, region in enumerate(data):
        for j, p in enumerate(region['points']):
            v = np.array(list(p))
            region['points'][j] = tuple(rotation.dot(v).tolist())
        region['center'] = rotation.dot(region['center'])

    # Rotate the coordinates around the z-axis 45 degrees
    #
    # R =  [ cos(45) -sin(45) ] = [ 1/sqrt(2)    -1/sqrt(2) ]
    #      [ sin(45)  sin(45) ]   [ 1/sqrt(2)   1/sqrt(2) ]
    rotation = (1 / np.sqrt(2)) * np.array([
        [1, -1, 0],
        [1, 1, 0],
        [0, 0, np.sqrt(2)]
    ])

    for i, region in enumerate(data):
        for j, p in enumerate(region['points']):
            v = np.array(list(p))
            region['points'][j] = tuple(rotation.dot(v).tolist())
        region['center'] = rotation.dot(region['center'])

    # Re-center
    device_center = np.array([x_center, y_center, 0])
    for i, region in enumerate(data):
        for j, p in enumerate(region['points']):
            v = np.array(list(p))
            region['points'][j] = tuple((v + device_center).tolist())
        region['center'] = device_center + region['center']

    return data


def get_triangle_regions(shunt_depth: float = 1, n_regions: int = 100, x_center: float = 0,
                         y_center: float = 0):
    """
    Defines the shunt regions for a 3D device

    Parameters
    ----------
    shunt_depth: float
        The depth of the shunt in um
    n_regions: int
        The number of conductivity regions
    x_center: float
        The x position of the base side of the shunt (um)
    y_center: float
        The x position of the base side of the shunt (um)

    Returns
    -------

    """
    # The side of the triangle
    L = np.sqrt(2.) * shunt_depth
    dy = np.sqrt(3.) * L / n_regions / 2.
    dx = L / 2 / n_regions
    z_offset = 0

    data = []
    x_l, x_r = - 0.5 * L, 0.5 * L
    isqrt3 = 1.0 / np.sqrt(3.)
    sqrt_23 = np.sqrt(2. / 3.)
    y = 0
    point_dtype = np.dtype([('x', 'd'), ('y', 'd'), ('z', 'd')])
    for i in range(n_regions):
        x_l_next, x_r_next, y_next = x_l + dx, x_r - dx, y + dy
        if i < n_regions - 1:
            points = np.empty(5, dtype=point_dtype)
            points[0] = (x_l, y, z_offset)
            points[1] = (x_r, y, z_offset)
            points[2] = (x_r_next, y_next, z_offset)
            points[3] = (x_l_next, y_next, z_offset)
            points[4] = (x_l, y, z_offset)
            data.append({
                'region': i,
                'points': points,
                'center': np.array([0.5 * (x_r + x_l), 0.5 * (y_next + y), z_offset])
            })
        else:
            points = np.empty(4, dtype=point_dtype)
            points[0] = (x_l, y, z_offset)
            points[1] = (x_r, y, z_offset)
            points[2] = (0, np.sqrt(3.) * L / 2., z_offset)
            points[3] = (x_l, y, z_offset)
            data.append({
                'region': i,
                'points': points,
                'center': np.array([0.5 * (x_r + x_l), 0.5 * (y_next + y), z_offset])
            })
        x_l, x_r, y = x_l_next, x_r_next, y_next

    # Rotate the coordinates around the x-axis 54 degrees
    #      
    # R =  [ cos(-54.73) -sin(-54.73) ] = [ 1/sqrt(3)    sqrt(2/3) ]
    #      [ sin(-54.73)  sin(-54.73) ]   [ -sqrt(2/3)   1/sqrt(3) ]
    rotation = np.array([
        [1., 0., 0],
        [0., isqrt3, sqrt_23],
        [0., -sqrt_23, isqrt3]
    ])
    for i, region in enumerate(data):
        for j, p in enumerate(region['points']):
            v = np.array(list(p))
            region['points'][j] = tuple(rotation.dot(v).tolist())
        region['center'] = rotation.dot(region['center'])

    # Rotate the coordinates around the z-axis 45 degrees
    #
    # R =  [ cos(45) -sin(45) ] = [ 1/sqrt(2)    -1/sqrt(2) ]
    #      [ sin(45)  sin(45) ]   [ 1/sqrt(2)   1/sqrt(2) ]
    rotation = (1 / np.sqrt(2)) * np.array([
        [1, -1, 0],
        [1, 1, 0],
        [0, 0, np.sqrt(2)]
    ])

    for i, region in enumerate(data):
        for j, p in enumerate(region['points']):
            v = np.array(list(p))
            region['points'][j] = tuple(rotation.dot(v).tolist())
        region['center'] = rotation.dot(region['center'])

    # Re-center
    device_center = np.array([x_center, y_center, 0])
    for i, region in enumerate(data):
        for j, p in enumerate(region['points']):
            v = np.array(list(p))
            region['points'][j] = tuple((v + device_center).tolist())
        region['center'] = device_center + region['center']

    return data


def conductivity_model(concentration: np.ndarray, segregation_coefficient: float = 1,
                       activated_na: float = 1) -> np.ndarray:
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
        segregation_coefficient: float
            The value of the segregation coefficient (unitless)

        Returns
        -------
        np.ndarray
            The conductivity_model profile
        """

    # Na concentration in the shunt
    cshunt = concentration * segregation_coefficient * activated_na

    # Model based on implantation data
    # Korol, V. M. "Sodium ion implantation into silicon." Physica status solidi (a) 110.1 (1988): 9-34.
    # Fitting of coefficients in Extract_NaImp.py
    coord = -11.144769029961262
    slope = 0.717839509854622

    sigma = np.power(10, coord) * np.power(cshunt, slope)  # (10 ** coord) * (cshunt ** slope)  # S/cm

    return sigma
