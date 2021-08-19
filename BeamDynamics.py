import numpy as np



C = 2.998e8   # Speed of light in [m/s]
E_REST_ELECTRON = 0.5100   # Rest energy of the electron [MeV]


def z_to_t(z, pz, Erest):
    beta = np.sqrt(pz**2. / (pz**2. + (Erest/c)**2.))
    t = z / (beta*c)
    return t
