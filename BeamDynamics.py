import os
import numpy as np
import pandas as pd
import json



C = 2.998e8   # Speed of light in [m/s]
E_REST_ELECTRON = 0.5100   # Rest energy of the electron [MeV]

COLUMN_ORDER_STANDARD_DF = ['x', 'px', 'y', 'py', 'z', 'pz', 't', 'E', 'gammaRel', 'betaRel', 'xp', 'yp']
UNITS_STANDARD_DF = ['mm', 'MeV/x', 'mm', 'MeV/c', 'mm', 'MeV/c', 'ns', 'MeV', '', '', 'mrad', 'mrad']


def p_to_beta(p, Erest):
    beta = np.sqrt(p**2. / (p**2. + (Erest/C)**2.))
    return beta


def z_to_t(z, pz, Erest):
    beta = p_to_beta(pz, Erest)
    t = z / (beta*C)
    return t


def pVect_to_p(px, py, pz):
    p = np.sqrt(px**2. + py**2. + pz**2.)
    return p


def pTransv_to_slope(pTransv, pLong):
    slope = pTransv / pLong * 1e3 # [mrad]
    return slope


def p_to_E(p, Erest):
    E = np.sqrt(p**2. + Erest**2.)
    return E


def E_to_gamma(E, Erest):
    gamma = E / Erest
    return gamma


def export_sim_db_to_xlsx(jsonFilePath):
    dbFile = open(jsonFilePath)
    dbJson = json.load(dbFile)
    dbDf = pd.json_normalize(dbJson)
    xlsxFilePath = os.path.splitext(jsonFilePath)[0] + '.xlsx'
    dbDf.to_excel(xlsxFilePath)


def convert_irina_distr_to_standard_df(sourceFilePath, saveStandardCsv=False):
    standardDf = pd.read_csv(
        sourceFilePath, delim_whitespace=True, index_col=False,
        header=0, names=('x', 'px', 'y', 'py', 'pz', 't')
    )
    # TODO: Implement exact conversion taking the relativistic beta into account
    standardDf['z'] = standardDf['t']
    # TODO: Choose unit for t
    standardDf['t'] = standardDf['t'] / C * 1e6 # [ns]
    p = pVect_to_p(standardDf['px'], standardDf['py'], standardDf['pz'])
    standardDf['E'] = p_to_E(p, E_REST_ELECTRON)
    standardDf['gammaRel'] = E_to_gamma(standardDf['E'], E_REST_ELECTRON)
    standardDf['betaRel'] = p_to_beta(standardDf['pz'], E_REST_ELECTRON)
    standardDf['xp'] = pTransv_to_slope(standardDf['px'], standardDf['pz'])
    standardDf['yp'] = pTransv_to_slope(standardDf['py'], standardDf['pz'])
    standardDf = standardDf[COLUMN_ORDER_STANDARD_DF]
    if saveStandardCsv:
        # headerList = standardDf.columns + ' [' + UNITS_STANDARD_DF + ']'
        standardTxtFilePath = os.path.splitext(sourceFilePath)[0] + '_StandardDf.txt'
        standardDf.to_string(standardTxtFilePath, index=False)
        # standardDf.to_csv(standardCsvFilePath, header=headerList)
    return standardDf