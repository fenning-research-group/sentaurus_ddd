import os
from typing import List


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
    parent_path = os.path.dirname(path)
    if not os.path.exists(parent_path):
        raise FileNotFoundError('Could not find path: \'{0}\'.'.format(parent_path))
    basename = os.path.basename(path)
    folder_name = next_file_name(parent_path, basename)
    full_path = os.path.join(parent_path, folder_name)
    os.makedirs(full_path)
    return full_path
