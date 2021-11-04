import numpy as np


                            
def compute_pcubed_bending_angle(zB1, zB2, alphaBuilding, dxCorner, zCorner, wBunker):
    A = zB2 - zB1
    B = zCorner + np.sin(np.deg2rad(alphaBuilding)) * wBunker/2. - zB2
    C = dxCorner - np.cos(np.deg2rad(alphaBuilding)) * wBunker/2.
    t = np.roots([A, -C, -A-2.*B, C])
    alphaB1 = np.rad2deg(np.arctan(t))
    return alphaB1

