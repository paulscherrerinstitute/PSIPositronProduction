import os
import subprocess
import warnings
import numpy as np
import pandas as pd
try:
    import ROOT
except:
    print('ROOT module not available.')
import json
import matplotlib.pyplot as plt



C = 2.99792458e8   # Speed of light in [m/s]
PART_CONSTS = {   # Particle constants
    'Erest': {   # Rest energy [MeV]
        11: 0.510998950,   # electron
        -11: 0.510998950,   # positron
        22: 0   # photon
    },
    'Q': {   # Electric charge [C]
        11: -1.60217663e-19,   # electron
        -11: 1.60217663e-19,   # positron
        22: 0   # photon
    }
}

COLUMN_ORDER_STANDARD_DF = ['x', 'px', 'y', 'py', 'z', 'pz', 't', 'E', 'gammaRel', 'betaRel', 'xp', 'yp', 'pdgId', 'Q', 'trackingId']
UNITS_STANDARD_DF = {
    'x': 'mm', 'px': 'MeV/c',
    'y': 'mm', 'py': 'MeV/c',
    'z': 'mm', 'pz': 'MeV/c',
    't': 'ps', 'E': 'MeV',
    'gammaRel': '', 'betaRel': '',
    'xp': 'mrad', 'yp': 'mrad',
    'pdgId': '', 'Q': 'C',
    'trackingId': ''
}
PRECISION_STANDARD_DF = 9
#PRECISION_STANDARD_DF = [6, 6, 6, 6, 6, 6, 9, 6, 3, 9, 6, 6]
#COLUMN_WIDTH_STANDARD_DF = 16

DATA_BASE_PATH = '/afs/psi.ch/project/Pcubed'


def build_data_path(relPath):
    if relPath[0] in ['/', '\\']:
        raise ValueError('relPath cannot start with a path separator ({:s} found).'.format(relPath[0]))
    absPath = os.path.join(DATA_BASE_PATH, relPath)
    return absPath


def pdgId_to_particle_const(pdgId, constName):
    Erest = pd.Series([
        PART_CONSTS[constName][pid] if pid in PART_CONSTS[constName].keys() else np.nan for pid in pdgId
    ])
    return Erest


def p_to_beta(p, pdgId):
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    beta = np.sqrt(p**2. / (p**2. + Erest**2.))
    return beta


def z_to_t(z, pz, pdgId):
    beta = p_to_beta(pz, pdgId)
    t = z / (beta*C)
    return t


def pVect_to_p(px, py, pz):
    p = np.sqrt(px**2. + py**2. + pz**2.)
    return p


def pTransv_to_slope(pTransv, pLong):
    slope = pTransv / pLong * 1e3 # [mrad]
    return slope


def p_to_E(p, pdgId):
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    E = np.sqrt(p**2. + Erest**2.)
    return E


def E_to_gamma(E, pdgId):
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    gamma = E / Erest
    return gamma


def export_sim_db_to_xlsx(jsonFilePath):
    dbFile = open(jsonFilePath)
    dbJson = json.load(dbFile)
    dbDf = pd.json_normalize(dbJson)
    xlsxFilePath = os.path.splitext(jsonFilePath)[0] + '.xlsx'
    dbDf.to_excel(xlsxFilePath)


def save_standard_fwf(sourceFilePath, standardDf):
    standardTxtFilePath = sourceFilePath + '.sdf_txt'
    formatterStr = '{:'+ str(PRECISION_STANDARD_DF+9) + '.' + str(PRECISION_STANDARD_DF) + 'e}'
    formatters = [formatterStr.format] * len(COLUMN_ORDER_STANDARD_DF)
    #formatters={
    #    'x': '',
    #    'pz': '{:'+str(COLUMN_WIDTH_STANDARD_DF)+'.6f}'.format,
    #}
    standardDf.to_string(standardTxtFilePath, index=False, formatters=formatters)
    # headerList = standardDf.columns + ' [' + UNITS_STANDARD_DF + ']'
    # standardDf.to_csv(standardCsvFilePath, header=headerList)


def load_standard_fwf(sourceFilePath):
    standardDf = pd.read_fwf(sourceFilePath)
    return standardDf


def convert_irina_distr_to_standard_df(sourceFilePath, saveStandardFwf=False):
    standardDf = pd.read_csv(
        sourceFilePath, delim_whitespace=True, index_col=False,
        header=0, names=('x', 'px', 'y', 'py', 'pz', 't')
    )
    standardDf['z'] = 17.5   # [mm]
    standardDf['t'] = standardDf['t'] / C * 1e9
    standardDf['pdgId'] = -11
    p = pVect_to_p(standardDf['px'], standardDf['py'], standardDf['pz'])
    standardDf['E'] = p_to_E(p, standardDf['pdgId'])
    standardDf['gammaRel'] = E_to_gamma(standardDf['E'], standardDf['pdgId'])
    standardDf['betaRel'] = p_to_beta(p, standardDf['pdgId'])
    standardDf['xp'] = pTransv_to_slope(standardDf['px'], standardDf['pz'])
    standardDf['yp'] = pTransv_to_slope(standardDf['py'], standardDf['pz'])
    standardDf['Q'] = pdgId_to_particle_const(standardDf['pdgId'], 'Q')
    standardDf = standardDf[COLUMN_ORDER_STANDARD_DF]
    if saveStandardFwf:
        save_standard_fwf(sourceFilePath, standardDf)
    return standardDf


def convert_fcceett_to_standard_df(sourceFilePath, pdgId=[-11], saveStandardFwf=False):
    if not isinstance(pdgId, list):
        pdgId = [pdgId]
    rootFile = ROOT.TFile.Open(sourceFilePath, 'READ')
    dfDict = {}
    for rootKey in rootFile.GetListOfKeys():
        distrName = rootKey.GetName()
        rootTree = rootFile.Get(distrName)
        rootTreeRdf = ROOT.RDataFrame(rootTree)
        rootTreeDict = rootTreeRdf.AsNumpy()
        rootTreeMat = np.column_stack([q for q in rootTreeDict.values()])
        standardDf = pd.DataFrame(data=rootTreeMat, columns=rootTreeDict.keys())
        standardDf.drop('evtId', axis=1, inplace=True)
        standardDf.rename(columns={'e': 'E'}, inplace=True)
        standardDf['pdgId'] = standardDf['pdgId'].astype(int)
        standardDf['gammaRel'] = E_to_gamma(standardDf['E'], standardDf['pdgId'])
        p = pVect_to_p(standardDf['px'], standardDf['py'], standardDf['pz'])
        standardDf['betaRel'] = p_to_beta(p, standardDf['pdgId'])
        standardDf['xp'] = pTransv_to_slope(standardDf['px'], standardDf['pz'])
        standardDf['yp'] = pTransv_to_slope(standardDf['py'], standardDf['pz'])
        standardDf['Q'] = pdgId_to_particle_const(standardDf['pdgId'], 'Q')
        standardDf = standardDf[COLUMN_ORDER_STANDARD_DF]
        fileSuffix = '_' + distrName
        if distrName == 'amor_leave' and pdgId:
            standardDf = standardDf[standardDf['pdgId'].isin(pdgId)]
            fileSuffix += '_pdgId'
            for id in pdgId:
                fileSuffix += '_' + str(id)
        if saveStandardFwf:
            filePath, fileExt = os.path.splitext(sourceFilePath)
            save_standard_fwf(filePath+fileSuffix+fileExt, standardDf)
        dfDict[distrName] = standardDf
    return dfDict


def convert_sdds_to_standard_df(sourceFilePath, z=np.nan, pdgId=-11, Qbunch=np.nan, saveStandardFwf=False):
    os.system('sddsconvert -ascii ' + sourceFilePath)
    standardDf = pd.read_csv(
        sourceFilePath, skiprows=24, delim_whitespace=True,
        names=('x', 'xp', 'y', 'yp', 't', 'p', 'trackingId')
    )
    pCentral = get_sdds_parameter(sourceFilePath, 'pCentral')
    QbunchRead = get_sdds_parameter(sourceFilePath, 'Charge')
    if not np.isnan(Qbunch) and not np.isclose(Qbunch, QbunchRead):
        warnStr = "Declared Qbunch = {:.6g} C (function input) does not correspond to value in file QbunchRead = {:.6f} C. Going to use ".format(
                Qbunch, QbunchRead
        )
        if not np.isnan(QbunchRead):
            Qbunch = QbunchRead
            warnStr += 'QbunchRead.'
        else:
            warnStr += 'declared Qbunch (function input).'
        warnings.warn(warnStr)
    else:
        Qbunch = QbunchRead
    NparticlesPerBunch = get_sdds_parameter(sourceFilePath, 'Particles')
    standardDf['x'] = standardDf['x'] * 1.e3                                    # [mm]
    standardDf['px'] = standardDf['xp'] * pCentral                              # [MeV/c]
    standardDf['xp'] = standardDf['xp'] * 1.e3                                  # [mrad]
    standardDf['y'] = standardDf['y'] * 1.e3                                    # [mm]
    standardDf['py'] = standardDf['yp'] * pCentral                              # [MeV/c]
    standardDf['yp'] = standardDf['yp'] * 1.e3                                  # [mrad]
    standardDf['z'] = z                                                         # [mm]
    standardDf['pdgId'] = pdgId
    p = standardDf['p'] \
        * pdgId_to_particle_const(standardDf['pdgId'], 'Erest')                 # [MeV/c]
    standardDf.drop('p', axis=1, inplace=True)
    standardDf['E'] = p_to_E(p, standardDf['pdgId'])                            # [MeV]
    standardDf['gammaRel'] = E_to_gamma(standardDf['E'], standardDf['pdgId'])
    standardDf['betaRel'] = p_to_beta(p, standardDf['pdgId'])
    standardDf['pz'] = \
        np.sqrt(p**2. - standardDf['px']**2. - standardDf['py']**2.)            # [MeV/c]
    standardDf['t'] = standardDf['t'] * 1.e9                                    # [ps]
    standardDf['Q'] = Qbunch / NparticlesPerBunch                               # [C]
    standardDf = standardDf[COLUMN_ORDER_STANDARD_DF]
    if saveStandardFwf:
        save_standard_fwf(sourceFilePath, standardDf)
    return standardDf


def get_sdds_parameter(sourceFilePath, parameterName):
    result = subprocess.run(
        ['sdds2stream', sourceFilePath, '-parameters='+parameterName],
        stdout=subprocess.PIPE
    )
    resultStr = result.stdout.decode('utf-8')
    try:
        resultVal = float(resultStr)
    except ValueError:
        resultVal = np.nan
    if np.isnan(resultVal):
        warnings.warn(
            "Could not read parameter {:s} correctly (resultStr = '{:s}').".format(
            parameterName, resultStr
        ))
    return resultVal


def plot_hist(ax, distr, binWidth, binLims=None, orientation='vertical', parsInLabel=True, opacity=1.):
    if binLims is None:
        binLims = [np.min(distr), np.max(distr) + binWidth]
    binEdges = np.arange(binLims[0], binLims[1], binWidth)
    counts, _, histObj = ax.hist(distr, bins=binEdges, orientation=orientation, alpha=opacity)
    std = np.std(distr)
    if parsInLabel:
        histObj.set_label('std = {:.3f}'.format(std))
        ax.legend()
    return std


def set_lims(ax, direction, var, lims):
    if lims is None:
        lims = (np.min(var), np.max(var))
    if not np.isclose(*lims):
        if direction == 'x':
            ax.set_xlim(lims)
        elif direction == 'y':
            ax.set_ylim(lims)


def plot_phase_space_2d(ax, distr, varName1=None, varName2=None, title=None, binWidth1=None, binWidth2=None, lims1=None, lims2=None, pzCutoff=None, opacity=1.):
    if varName1 is None or varName2 is None:
        raise ValueError('varName1 and varName2 are mandatory arguments.')
    markerSize = 10
    ax[0,0].scatter(distr[varName1], distr[varName2], markerSize, alpha=opacity)
    plot_hist(ax[1,0], distr[varName1], binWidth1, binLims=lims1, opacity=opacity)
    plot_hist(ax[0,1], distr[varName2], binWidth2, binLims=lims2, orientation='horizontal', opacity=opacity)
    pzLabels = ['All, Counts = {:d}'.format(distr.shape[0]), ]
    if pzCutoff is not None:
        distrLowPz = distr[distr['pz']<pzCutoff]
        ax[0,0].scatter(distrLowPz[varName1], distrLowPz[varName2], markerSize, alpha=opacity)
        plot_hist(ax[1,0], distrLowPz[varName1], binWidth1, binLims=lims1, opacity=opacity)
        plot_hist(ax[0,1], distrLowPz[varName2], binWidth2, binLims=lims2, orientation='horizontal', opacity=opacity)
        pzLabels += ['pz <= {:.1f} {:s}, Counts = {:d}'.format(
            pzCutoff, UNITS_STANDARD_DF['pz'], distrLowPz.shape[0]
        )]
        ax[0,0].legend(pzLabels, markerscale=5., loc='upper left', bbox_to_anchor=(1.2, -0.2))
    set_lims(ax[0,0], 'x', distr[varName1], lims1)
    set_lims(ax[0,0], 'y', distr[varName2], lims2)
    ax[1,0].set_xlim(ax[0,0].get_xlim())
    ax[1,0].invert_yaxis()
    ax[0,1].set_ylim(ax[0,0].get_ylim())
    ax[0,0].grid()
    ax[0,1].grid()
    ax[1,0].grid()
    ax[0,0].set_xlabel(varName1+' ['+UNITS_STANDARD_DF[varName1]+']')
    ax[0,0].set_ylabel(varName2+' ['+UNITS_STANDARD_DF[varName2]+']')
    ax[0,0].yaxis.set_label_position('right')
    ax[0,0].yaxis.tick_right()
    ax[1,0].set_ylabel('Bin counts')
    ax[0,1].set_xlabel('Bin counts')
    if title is not None:
        ax[0,0].get_figure().suptitle(title)
    ax[1,1].set_visible(False)


def plot_distr(distr, plotDefs, title=None, figHeight=9, figWidth=16):
    for d in plotDefs:
        fig, ax = plt.subplots(2, 2)
        fig.set_figheight(figHeight)
        fig.set_figwidth(figWidth)
        plot_phase_space_2d(ax, distr, **d, title=title)
        plt.show()
