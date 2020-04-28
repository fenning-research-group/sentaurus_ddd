import numpy as np
import h5py
import pandas as pd
import os
from typing import List

region_attr_type = np.dtype(
    [
        ('region_id', 'U100'),
        ('region index', 'i'),
        ('material', 'U100'),
        ('material type', 'U100'),
        ('name', 'U100'),
        ('number of parts', 'i'),
        ('type', 'i')
    ]
)

# Structure type: 0 = scalar, 1 = vector: 2 = tensor
# location type: 0 = vertex centered, 1 = element centered
# value type: 0 = integer, 1 = float, 2 = double
dataset_attr_type = np.dtype(
    [
        ('dataset_id', 'U100'),
        ('dataset index', 'i'),
        ('conversion factor', 'd'),
        ('location type', 'i'),
        ('name', 'U100'),
        ('number of values', 'i'),
        ('quantity', 'U100'),
        ('region index', 'i'),
        ('structure type', 'i'),
        ('unit:long name', 'U100'),
        ('unit:name', 'U100'),
        ('unit:significand', 'd'),
        ('value type', 'i')
    ]
)


class TDR:
    """
    This class provides access to 2D (and 3D) data sets in Sentaurus TDR files
    """

    _dataset_tree: pd.DataFrame = None
    _region_tree: pd.DataFrame = None
    _space_dimension: int = 2
    _vertices_per_element: int = 3

    def __init__(self, tdr_file: str, space_dimension: int = 2):
        self._tdr_file = tdr_file
        self._vertex = self.get_vertex()
        self._region_tree = self.get_region_tree()
        self._dataset_tree = self.get_dataset_tree()
        self._space_dimension = space_dimension
        if space_dimension == 2:
            self._vertices_per_element = 3
        elif space_dimension == 3:
            self._vertices_per_element = 4

    def get_dataset_ids(self, name: str, region_idx: str = None, material: str = None):
        datasets_df = self._dataset_tree[self._dataset_tree['name'] == name]
        if material is not None:
            regions_df = self.get_region_by_material(material)
            return pd.merge(datasets_df, regions_df, on='region index', how='inner').reset_index(drop=True)
        elif region_idx is not None:
            regions_df = self._region_tree[self._region_tree['region index'] == region_idx]
            return pd.merge(datasets_df, regions_df, on='region index', how='inner').reset_index(drop=True)
        return datasets_df

    def get_dataset(self, name: str, region_idx: str = None, material: str = None):
        ds_info = self.get_dataset_ids(name=name, region_idx=region_idx, material=material)
        # Check whether the dataset is vertex centered or element centered
        location_type_int = ds_info.iloc[0]['location type']
        location_type = 'vertex-centered' if location_type_int == 0 else 'element-centered'
        # If location is vertex centered
        if location_type == 'vertex-centered':
            # Get all the vertex
            vertex_all = self.get_vertex()
            # Get all the vertices in the dataset
            vertex_indices = []
            values = []
            for i, r in ds_info.iterrows():
                region_idx: int = r['region index']
                vertex_indices = self.get_region_vertex_indices(region_idx=region_idx)
                for v in vertex_indices:
                    coords = vertex_all
                    vertex_indices.append()

    def get_region_vertex_indices(self, region_idx: int):
        group_name = 'collection/geometry_0/region_{0:d}/elements_0'.format(region_idx)
        with h5py.File(self._tdr_file, 'r') as hf:
            elements = np.array(hf.get(group_name))
            # If dimension = 2D each element should look like (type, index1, index2, index3)
            if self._space_dimension == 2:
                vertex_per_element = 3
                n_elements = int(len(elements) / 4)
            # If dimension = 3D each element should look like (type, index1, index2, index3, index4)
            elif self._space_dimension == 3:
                vertex_per_element = 4
                n_elements = int(len(elements) / 5)
            elements = elements.reshape((n_elements, vertex_per_element))
            vertex_ids = elements[:,1::]
            unique_vertex = []
            for v in vertex_ids:
                if v not in unique_vertex:
                    unique_vertex.append(v)
            return unique_vertex

    def get_vertex(self):
        group_name = 'collection/geometry_0/vertex'
        with h5py.File(self._tdr_file, 'r') as hf:
            vertex_all = np.array(hf[group_name])
        return vertex_all

    def get_region_by_material(self, material: str) -> pd.DataFrame:
        return self.get_region_by_attribute(attribute_name='material', value=material)

    def get_region_by_attribute(self, attribute_name: str, value) -> pd.DataFrame:
        return self._region_tree[self._region_tree[attribute_name] == value].reset_index(drop=True)

    def get_region_tree(self) -> pd.DataFrame:
        group_name = 'collection/geometry_0'
        with h5py.File(self._tdr_file, 'r') as hf:
            group = hf.get(group_name)
            regions = [r for r in list(group) if r.startswidth('region')]
            n_regions = len(regions)
            region_tree = np.empty(n_regions, dtype=region_attr_type)
            for i, r in enumerate(regions):
                region_attrs = self.region_attributes(region_name=r)
                region_tree[i] = (r, i, region_attrs['material'].decode('utf-8'),
                                  region_attrs['material type'].decode('utf-8'),
                                  region_attrs['name'].decode('utf-8'),
                                  int(region_attrs['number of parts']),
                                  int(region_attrs['type']))
        return pd.DataFrame(data=region_tree)

    def get_dataset_tree(self) -> pd.DataFrame:
        group_name = 'collection/geometry_0/state_0'
        with h5py.File(self._tdr_file, 'r') as hf:
            group = hf.get(group_name)
            datasets = [r for r in list(group) if r.startswidth('dataset')]
            n_datasets = len(datasets)
            dataset_tree = np.empty(n_datasets, dtype=dataset_attr_type)
            for i, r in enumerate(datasets):
                dataset_attrs = self.dataset_attributes(dataset_name=r)
                dataset_attrs[i] = (r, i, dataset_attrs['conversion factor'],
                                    dataset_attrs['location type'],
                                    dataset_attrs['number of values'],
                                    dataset_attrs['quantity'].decode('utf-8'),
                                    dataset_attrs['region'],
                                    dataset_attrs['structure type'],
                                    dataset_attrs['unit:long name'].decode('utf-8'),
                                    dataset_attrs['unit:name'].decode('utf-8'),
                                    dataset_attrs['unit:significand'],
                                    dataset_attrs['value type'])
        return pd.DataFrame(data=dataset_tree)

    def region_attributes(self, region_name: str):
        return self.get_group_attributes(group_name='collection/geometry_0/{0}'.format(region_name))

    def dataset_attributes(self, dataset_name):
        return self.get_group_attributes(group_name='collection/geometry_0/state_0/{0}'.format(dataset_name))

    def dataset_by_name(self, name: str):
        with h5py.File(self._tdr_file, 'r') as hf:
            group = hf[name]
        return np.array(group)

    def get_attribute(self, name: str, group: str = None):
        with h5py.File(self._tdr_file, 'r') as hf:
            if group is None:
                value = hf.attrs[name]
            else:
                g = hf.get(group)
                value = g.attrs[name]
        return value

    def get_group_attributes(self, group_name: str) -> dict:
        with h5py.File(self._tdr_file, 'r') as hf:
            g = hf[group_name]
            attributes = {a: g.attrs[a] for a in g.attrs}
        return attributes

    def get_vertex(self):
        with h5py.File(self._tdr_file, 'r') as hf:
            v = np.array(hf['collection/geometry_0/vertex'])
        return v
