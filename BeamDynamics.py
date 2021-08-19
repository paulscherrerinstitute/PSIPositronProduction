import os
import numpy as np
import pandas as pd
import json



C = 2.998e8   # Speed of light in [m/s]
E_REST_ELECTRON = 0.5100   # Rest energy of the electron [MeV]


def z_to_t(z, pz, Erest):
    beta = np.sqrt(pz**2. / (pz**2. + (Erest/C)**2.))
    t = z / (beta*C)
    return t


def export_sim_db_to_xlsx(jsonFilePath):
    dbFile = open(jsonFilePath)
    dbJson = json.load(dbFile)
    dbDf = pd.json_normalize(dbJson)
    xlsxFilePath = os.path.splitext(jsonFilePath)[0] + '.xlsx'
    dbDf.to_excel(xlsxFilePath)
