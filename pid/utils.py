import os
from typing import List
import numpy as np
import pandas as pd


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


def geometric_series_spaced(max_val: float, min_delta: float, steps: int, **kwargs) -> np.ndarray:
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
    **kwargs: keyword arguments
        n_iterations: int
            The number of Newton iterations to estimate r

    Returns
    -------
    np.ndarray
        A vector with geometrically spaced values
    """
    niter = kwargs.get('n_iterations', 1000)
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
    # Estimate r
    r = 2
    for i in range(niter):
        if fp(r) == 0:  # Avoid zeros
            r += 0.1
        r -= f(r)/fp(r)

    result = np.zeros(n + 1)
    s = 0
    p = 1
    for i in range(n):
        s += a*p
        result[i+1] = s
        p *= r

    return result


def get_indices_at_values(x: np.array, requested_values: np.array) -> List[int]:
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
    L = np.sqrt(2.)*shunt_depth
    dy = np.sqrt(3.)*L/n_regions/2.
    dx = L / 2 / n_regions
    z_offset = 0

    data = []
    x_l, x_r = - 0.5 * L,  0.5 * L
    isqrt3 = 1.0/np.sqrt(3.)
    sqrt_23 = np.sqrt(2./3.)
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
                'center': np.array([0.5*(x_r + x_l), 0.5*(y_next+y), z_offset])
            })
        else:
            points = np.empty(4, dtype=point_dtype)
            points[0] = (x_l, y, z_offset)
            points[1] = (x_r, y, z_offset)
            points[2] = (0, np.sqrt(3.)*L/2., z_offset)
            points[3] = (x_l, y, z_offset)
            data.append({
                'region': i,
                'points': points,
                'center': np.array([0.5 * (x_r + x_l), 0.5*(y_next+y), z_offset])
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
            # print(region['points'][j])
        # p1 = np.array(region['point1']).T
        # p2 = np.array(region['point2']).T
        # p3 = np.array(region['point3']).T
        #
        # region['point1'] = rotation.dot(p1)
        # region['point2'] = rotation.dot(p2)
        # region['point3'] = rotation.dot(p3)
        #
        # if i < n_regions - 1:
        #     p4 = np.array(region['point4']).T
        #     region['point4'] = rotation.dot(p4)

    # Rotate the coordinates around the z-axis 45 degrees
    #
    # R =  [ cos(45) -sin(45) ] = [ 1/sqrt(2)    -1/sqrt(2) ]
    #      [ sin(45)  sin(45) ]   [ 1/sqrt(2)   1/sqrt(2) ]
    rotation = (1/np.sqrt(2))*np.array([
        [1, -1, 0],
        [1, 1, 0],
        [0, 0, np.sqrt(2)]
    ])

    for i, region in enumerate(data):
        for j, p in enumerate(region['points']):
            v = np.array(list(p))
            region['points'][j] = tuple(rotation.dot(v).tolist())
        region['center'] = rotation.dot(region['center'])
        # p1 = np.array(region['point1']).T
        # p2 = np.array(region['point2']).T
        # p3 = np.array(region['point3']).T
        #
        # region['point1'] = rotation.dot(p1)
        # region['point2'] = rotation.dot(p2)
        # region['point3'] = rotation.dot(p3)
        #
        # if i < n_regions - 1:
        #     p4 = np.array(region['point4']).T
        #     region['point4'] = rotation.dot(p4)

    # Re-center
    device_center = np.array([x_center, y_center, 0])
    for i, region in enumerate(data):
        for j, p in enumerate(region['points']):
            v = np.array(list(p))
            region['points'][j] = tuple((v + device_center).tolist())
        region['center'] = device_center + region['center']
        # p1 = np.array(region['point1']).T
        # p2 = np.array(region['point2']).T
        # p3 = np.array(region['point3']).T
        #
        # region['point1'] = p1 + device_center
        # region['point2'] = p2 + device_center
        # region['point3'] = p3 + device_center
        #
        # if i < n_regions - 1:
        #     p4 = np.array(region['point4']).T
        #     region['point4'] = p4 + device_center

    return data



