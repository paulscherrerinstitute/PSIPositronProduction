import numpy as np
import pandas as pd
import time

import os, sys
BIN = os.path.expanduser('~/Git/rf-track-2.1')
sys.path.append(BIN)

import warnings
warnings.filterwarnings('ignore')

def prepare_fieldmap(csv_filename, skiprows=8, convert_units=True, z_cut_m=None):
    "read a csv file and prepares for RF Track"

    fieldmap = pd.read_csv(csv_filename, delim_whitespace=True, skiprows=skiprows, names=['x', 'y', 'z', 'Bx', 'By', 'Bz'])

    # append electric field
    fieldmap[['Ex', 'Ey', 'Ez']] = pd.DataFrame([[0, 0, 0]], index = fieldmap.index)

    if convert_units:
        fieldmap[['x', 'y', 'z']] = fieldmap[['x', 'y', 'z']]*1e-3 # convert mm to m

    # remove part of the field map where abs(z) < z_cut_m
    if z_cut_m is None:
        z_cut_m = max(fieldmap.z)

    # indices of the field map dataframe to be removed
    assert z_cut_m >= 0, "Cut in z [m] must be greater than 0"
    indices_rm = fieldmap.index[np.abs(fieldmap['z']) > z_cut_m ].to_list()

    # field map with the boundaries
    fieldmap.loc[indices_rm, ['Bx', 'By', 'Bz']] = None, None, None
    fieldmap_bounds = fieldmap.dropna()

    # mesh step size
    mesh_stepsize = {}
    mesh_stepno = {}
    mesh_dimensions = {}
    mesh_bounds = {}
    for axis in ['x', 'y', 'z']:
        axis_val = fieldmap_bounds[axis].unique()

        mesh_stepno[axis] = len(axis_val)
        mesh_bounds[axis] = [min(fieldmap_bounds[axis]), max(fieldmap_bounds[axis])]
        mesh_dimensions[axis] = abs(min(fieldmap_bounds[axis]) - max(fieldmap_bounds[axis]))
        mesh_stepsize[axis] = np.unique(np.diff(axis_val)).mean().round(6)

    fieldmap_info = {'stepsizes': (mesh_stepsize['x'], mesh_stepsize['y'], mesh_stepsize['z']),
                     'bounds': (mesh_bounds['x'], mesh_bounds['y'], mesh_bounds['z']),
                     'dimensions': (mesh_dimensions['x'], mesh_dimensions['y'], mesh_dimensions['z']),
                     'stepnos': (mesh_stepno['x'], mesh_stepno['y'], mesh_stepno['z'])}

    return fieldmap_bounds, fieldmap_info

# Load RF-Track
import RF_Track as rft

# Load 3d field map of an example dipole (Mikko's dipole)
fieldmap_files = ["fieldmap_dx10dy10dz10_nx31ny11nz201.B",
                  "fieldmap_dx2dy2dz2_nx151ny51nz1001.B"]

fieldmap_file = fieldmap_files[1]

# read the field map
fieldmap, fieldmap_info = prepare_fieldmap(csv_filename=fieldmap_file, z_cut_m=0.7)

(nx, ny, nz) = fieldmap_info['stepnos']
(lx, ly, lz) = fieldmap_info['dimensions']
(hx, hy, hz) = fieldmap_info['stepsizes']
print(fieldmap_info['bounds'])

print("ciao = ", np.array(fieldmap.Bx).size)

Bx = np.array(fieldmap.Bx).reshape(nx, ny, nz)
By = np.array(fieldmap.By).reshape(nx, ny, nz)
Bz = np.array(fieldmap.Bz).reshape(nx, ny, nz)
fieldmap = 0
print("nx, ny, nz = ", nx, ny, nz)
full_dipole = rft.RF_FieldMap_CINT(0, 0, 0, \
                       Bx,  By, Bz, \
                       -lx/2, \
                       -ly/2, \
                       hx, \
                       hy, \
                       hz, \
                       lz, 0, 0)

full_dipole = 0

print("A")
