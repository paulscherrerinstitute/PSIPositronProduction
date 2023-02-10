import os
import subprocess
import warnings
import numpy as np
import pandas as pd
import scipy.stats as scistats
import scipy.interpolate as sciinterp
try:
    import ROOT
except ModuleNotFoundError:
    warnings.warn('ROOT module not available.')
import json
import matplotlib.pyplot as plt
import matplotlib.markers as pltMarkers
import matplotlib.patches as pltPatches
import OctavePythonInterface as opi


pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
# pd.set_option('display.max_colwidth', -1)


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

COLUMN_ORDER_STANDARD_DF = ['x', 'px', 'y', 'py', 'z', 'pz', 't', 'pdgId', 'Q']
UNITS_STANDARD_DF_EXTENDED = {
    'x': 'mm', 'px': 'MeV/c',
    'y': 'mm', 'py': 'MeV/c',
    'z': 'mm', 'pz': 'MeV/c',
    't': 'ns',
    'pdgId': '', 'Q': 'C',
    'Ekin': 'MeV', 'gammaRel': '', 'betaRel': '',
    'xp': 'mrad', 'yp': 'mrad',
    'r': 'mm'
}
PRECISION_STANDARD_DF = 9

SCI_NOTATION_FORMATTER = (
    '{:' + str(PRECISION_STANDARD_DF+9) + '.' + str(PRECISION_STANDARD_DF) + 'e}'
).format
INTEGER_FORMATTER = ('{:'+str(PRECISION_STANDARD_DF+9)+'d}').format

FILE_TYPES_SPECS = {
    'standardDf': {
        'ext': '.sdf_txt',
        'columnOrder': COLUMN_ORDER_STANDARD_DF,
        'header': True,
        'formatters': [SCI_NOTATION_FORMATTER]*len(COLUMN_ORDER_STANDARD_DF)
        # formatterStr = '{:'+ str(PRECISION_STANDARD_DF+9) + '.' + \
        # str(PRECISION_STANDARD_DF) + 'e}'
        # {formatters = [formatterStr.format] * fileTypeSpecs['numCols']
    },
    'astra': {
        'ext': '.001',
        'columnOrder': [
            'x', 'y', 'z', 'px', 'py', 'pz',
            'clock', 'macroCharge', 'particleIndex', 'statusFlag'
        ],
        'header': False,
        'formatters': {
            'x': SCI_NOTATION_FORMATTER,
            'y': SCI_NOTATION_FORMATTER,
            'z': SCI_NOTATION_FORMATTER,
            'px': SCI_NOTATION_FORMATTER,
            'py': SCI_NOTATION_FORMATTER,
            'pz': SCI_NOTATION_FORMATTER,
            'clock': SCI_NOTATION_FORMATTER,
            'macroCharge': SCI_NOTATION_FORMATTER,
            'particleIndex': INTEGER_FORMATTER,
            'statusFlag': INTEGER_FORMATTER
        }
    },
    'sdds': {
        'ext': '.sdds',
        'columnOrder': ['x', 'xp', 'y', 'yp', 't', 'p', 'trackingId'],
        'header': None
    },
    'rftrack_xp_t': {
        'ext': '.dat',
        'columnOrder': ['x', 'xp', 'y', 'yp', 't', 'p'],
        'columnUnits': ['mm', 'mrad', 'mm', 'mrad', 'mm/c', 'MeV/c'],
        'header': None,
        'formatters': {
            'x': SCI_NOTATION_FORMATTER,
            'xp': SCI_NOTATION_FORMATTER,
            'y': SCI_NOTATION_FORMATTER,
            'yp': SCI_NOTATION_FORMATTER,
            't': SCI_NOTATION_FORMATTER,
            'p': SCI_NOTATION_FORMATTER
        }
    },
    'rftrack_Px_S': {
        'ext': '.dat',
        # TODO: Is the last quantity really pz or p = sqrt(px^2+py^2+pz^2)?
        'columnOrder': ['x', 'px', 'y', 'py', 's', 'pz'],
        'columnUnits': ['mm', 'MeV/c', 'mm', 'MeV/c', 'mm', 'MeV/c'],
        'header': None,
        'formatters': {
            'x': SCI_NOTATION_FORMATTER,
            'px': SCI_NOTATION_FORMATTER,
            'y': SCI_NOTATION_FORMATTER,
            'py': SCI_NOTATION_FORMATTER,
            's': SCI_NOTATION_FORMATTER,
            'pz': SCI_NOTATION_FORMATTER
        }
    },
    'placet': {
        'ext': '.dat',
        'columnOrder': ['p', 'x', 'y', 't', 'xp', 'yp'],
        'header': False,
        'formatters': {
            'p': SCI_NOTATION_FORMATTER,
            'x': SCI_NOTATION_FORMATTER,
            'y': SCI_NOTATION_FORMATTER,
            't': SCI_NOTATION_FORMATTER,
            'xp': SCI_NOTATION_FORMATTER,
            'yp': SCI_NOTATION_FORMATTER,
        }
    },
}


def pdgId_to_particle_const(pdgId, constName):
    flagSingleValue = False
    if not isinstance(pdgId, (pd.Series, np.ndarray)):
        flagSingleValue = True
        pdgId = [pdgId]
    Erest = np.array([
        PART_CONSTS[constName][pid] if pid in PART_CONSTS[constName].keys() else np.nan
        for pid in pdgId
    ])
    if flagSingleValue:
        Erest = Erest[0]
    return Erest


def check_input_pd_or_np(inArray):
    try:
        outArray = inArray.to_numpy()
    except AttributeError:
        if isinstance(inArray, float):
            outArray = np.array(inArray)
        elif isinstance(inArray, np.ndarray):
            outArray = inArray
        else:
            raise TypeError("inArray must be either float, numpy.ndarray or pandas.Series.")
    return outArray


def p_to_beta(p, pdgId):
    p = check_input_pd_or_np(p)
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    beta = np.sqrt(p**2. / (p**2. + Erest**2.))
    return beta


def pVect_to_p(px, py, pz):
    px = check_input_pd_or_np(px)
    py = check_input_pd_or_np(py)
    pz = check_input_pd_or_np(pz)
    p = np.sqrt(px**2. + py**2. + pz**2.)
    return p


def pTransv_to_slope(pTransv, pLong):
    pTransv = check_input_pd_or_np(pTransv)
    pLong = check_input_pd_or_np(pLong)
    slope = np.arctan2(pTransv,  pLong)
    return slope * 1e3  # [mrad]


def p_to_Ekin(p, pdgId):
    p = check_input_pd_or_np(p)
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    E = np.sqrt(p**2. + Erest**2.)
    Ekin = E - Erest
    return Ekin


def Ekin_to_p(Ekin, pdgId):
    Ekin = check_input_pd_or_np(Ekin)
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    p = np.sqrt((Ekin+Erest)**2. - Erest**2.)
    return p


def Ekin_to_gamma(Ekin, pdgId):
    Ekin = check_input_pd_or_np(Ekin)
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    gamma = Ekin / Erest + 1
    return gamma


def gamma_to_beta(gamma):
    gamma = check_input_pd_or_np(gamma)
    beta = np.sqrt(1. - (1./gamma**2.))
    return beta


def gamma_to_p(gamma, pdgId):
    gamma = check_input_pd_or_np(gamma)
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    p = Erest * np.sqrt(gamma**2. - 1)
    return p


def p_to_gamma(p, pdgId):
    p = check_input_pd_or_np(p)
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    gamma = np.sqrt((p/Erest)**2. - 1.)
    return gamma


def p_to_betagamma(p, pdgId):
    p = check_input_pd_or_np(p)
    beta = p_to_beta(p, pdgId)
    gamma = p_to_gamma(p, pdgId)
    betaGamma = beta * gamma
    return betaGamma


def Ekin_to_beta(Ekin, pdgId):
    gamma = Ekin_to_gamma(Ekin, pdgId)
    beta = gamma_to_beta(gamma)
    return beta


def extend_standard_df(standardDf, removeNanInf):
    p = pVect_to_p(standardDf['px'], standardDf['py'], standardDf['pz'])
    dfExtension = pd.DataFrame(index=standardDf.index)
    dfExtension['Ekin'] = p_to_Ekin(p, standardDf['pdgId'])
    dfExtension['gammaRel'] = Ekin_to_gamma(dfExtension['Ekin'], standardDf['pdgId'])
    dfExtension['betaRel'] = p_to_beta(p, standardDf['pdgId'])
    dfExtension['xp'] = pTransv_to_slope(standardDf['px'], standardDf['pz'])
    dfExtension['yp'] = pTransv_to_slope(standardDf['py'], standardDf['pz'])
    # By definition, backwards traveling particles have xp > np.pi/2. or xp < -np.pi/2.
    # while -np.pi/2 < yp < np.pi/2 in standardDf
    # TODO: This definition should be implemented in every function generating a standardDf
    dfExtension.loc[(standardDf['pz'] < 0) & (dfExtension['yp'] > np.pi/2.*1e3), 'yp'] -= np.pi*1e3
    dfExtension.loc[(standardDf['pz'] < 0) & (dfExtension['yp'] < -np.pi/2.*1e3), 'yp'] += np.pi*1e3
    # Check conformity of input distribution to angle definition (if angles available)
    numDecimalPlaces = 6
    if 'yp' in standardDf.columns:
        indsAnglesNotConformal = (standardDf['yp'] > np.pi/2.*1e3) | \
            (standardDf['yp'] < -np.pi/2.*1e3)
        anglesNotConformal = standardDf.loc[indsAnglesNotConformal, ['xp', 'yp']]
        anglesNotConformal.loc[anglesNotConformal['yp'] > np.pi/2.*1e3, 'yp'] -= np.pi*1e3
        anglesNotConformal.loc[anglesNotConformal['yp'] < -np.pi/2.*1e3, 'yp'] += np.pi*1e3
        indsXpLargerZero = anglesNotConformal['xp'] > 0
        anglesNotConformal.loc[indsXpLargerZero, 'xp'] -= np.pi*1e3
        anglesNotConformal.loc[~indsXpLargerZero, 'xp'] += np.pi*1e3
        probAngles = dfExtension.loc[indsAnglesNotConformal, ['xp', 'yp']].round(numDecimalPlaces) \
            .compare(anglesNotConformal.round(numDecimalPlaces))
        if probAngles.shape[0] > 0:
            warnings.warn(
                '{:d} backward traveling particles '.format(probAngles.shape[0]) +
                'with problematic angles detected with index: ' +
                ', '.join('{:d}'.format(partInd) for partInd in probAngles.index)
            )
            print(probAngles.head())
        indsToCompare = standardDf.index[~indsAnglesNotConformal]
    else:
        indsToCompare = standardDf.index
    commonQuantities = standardDf.columns.intersection(dfExtension.columns)
    dfDiff = standardDf.loc[indsToCompare, commonQuantities].round(numDecimalPlaces) \
        .compare(dfExtension.loc[indsToCompare, commonQuantities].round(numDecimalPlaces))
    if dfDiff.shape[0] > 0:
        warnings.warn(
            '{:d} differing particles detected with index: '.format(dfDiff.shape[0]) +
            ', '.join('{:d}'.format(partInd) for partInd in dfDiff.index)
        )
        print(dfDiff)
    standardDf = dfExtension.combine_first(standardDf)
    standardDf = check_nan_inf_in_distr(standardDf, removeNanInf)
    return standardDf


def compute_emittance(
        standardDf, planeName, norm='normalized',
        correctOffsets=True, verbose=True,
        filterSpecs={}
):
    if norm in ['normalized', 'geometric']:
        uDivName = 'p' + planeName
        uDivUnits = 'MeV/c'
    elif norm == 'tracespace':
        uDivName = planeName + 'p'
        uDivUnits = 'mrad'
    else:
        raise ValueError(
            "norm must be either 'normalized', 'geometric' or 'tracespace'."
        )
    standardDf = filter_distr(standardDf, filterSpecs)
    # NOTE: Indexing a dataframe returns its reference, therefore use copy()!
    u = standardDf[planeName].copy()
    uDiv = standardDf[uDivName].copy()
    u, uDiv = check_distribution_offsets(
        u, uDiv, planeName, uDivName, uDivUnits, correctOffsets, verbose
    )
    # TODO: Check that all pdgIds are identical
    emit = np.sqrt(
        (u**2.).sum() * (uDiv**2.).sum() - (u*uDiv).sum()**2.
    ) / u.shape[0]
    if norm in ['normalized', 'geometric']:
        Erest = pdgId_to_particle_const(standardDf['pdgId'].iloc[0], 'Erest')
        emit *= 1e3 / Erest  # [mm mrad]
        if norm == 'geometric':
            emit /= (standardDf['betaRel']*standardDf['gammaRel']).mean()
    return emit  # [mm mrad]


def compute_twiss(
        standardDf, planeName, filterSpecs={},
        correctOffsets=True, verbose=True
):
    standardDf = filter_distr(standardDf, filterSpecs)
    emitTraceSpace = compute_emittance(
        standardDf, planeName, norm='tracespace',
        correctOffsets=correctOffsets, verbose=verbose
    )
    uDivName = planeName + 'p'
    uDivUnits = 'mrad'
    # TODO: refactor (see compute_emittance)
    # NOTE: Indexing a dataframe returns its reference, therefore use copy()!
    u, uDiv = check_distribution_offsets(
        standardDf[planeName].copy(), standardDf[uDivName].copy(),
        planeName, uDivName, uDivUnits, correctOffsets, verbose
    )
    alphaTwiss = -1. * (u*uDiv).sum()/u.shape[0] / emitTraceSpace
    betaTwiss = (u**2.).sum()/u.shape[0] / emitTraceSpace
    gammaTwiss = (uDiv**2.).sum()/u.shape[0] / emitTraceSpace
    return alphaTwiss, betaTwiss, gammaTwiss


def check_distribution_offsets(
        u, uDiv, planeName, uDivName, uDivUnits, correctOffsets, verbose
):
    thresholdFactor = 0.001
    if verbose or correctOffsets:
        uAvg = u.mean()
        uDivAvg = uDiv.mean()
    if verbose:
        if uAvg > u.std()*thresholdFactor:
            warnings.warn(
                'Average position {:s}Avg = {:.3f} mm.'.format(planeName, uAvg)
            )
        if uDivAvg > uDiv.std()*thresholdFactor:
            warnings.warn(
                'Average divergence {:s}Avg = {:.3f} {:s}.'.format(uDivName, uDivAvg, uDivUnits)
            )
    if correctOffsets:
        if not (np.isnan(uAvg) or np.isnan(uDivAvg)):
            u -= uAvg
            uDiv -= uDivAvg
            if verbose:
                print(
                    'Correcting offsets {:s}Avg = {:.3f} mm and {:s}Avg = {:.3f} {:s}.'
                    .format(planeName, uAvg, uDivName, uDivAvg, uDivUnits)
                )
        else:
            print(
                'Offset cannot be corrected with {:s}Avg = {:.3f} mm and {:s}Avg = {:.3f} {:s}.'
                .format(planeName, uAvg, uDivName, uDivAvg, uDivUnits)
            )
    return u, uDiv


def use_filter_specs_selector(standardDf, sourceFilePath, filterSpecsSelector):
    filterSpecs = get_json_entry(sourceFilePath, filterSpecsSelector, 'filterSpecs')
    standardDf = filter_distr(standardDf, filterSpecs)
    return standardDf


def filter_distr(standardDf, filterSpecs):
    for varName, varLims in filterSpecs.items():
        standardDf = standardDf[
            (standardDf[varName] >= varLims[0])
            & (standardDf[varName] <= varLims[1])
        ]
    return standardDf


def get_json_entry(sourceFilePath, specsSelector, propertyName):
    folderPath = os.path.dirname(sourceFilePath)
    # TODO: Rename standard file with better name?
    with open(os.path.join(folderPath, 'filterSpecs.json'), 'r') as specsFile:
        specsList = json.load(specsFile)
    selectedProperty = specsList[specsSelector][propertyName]
    return selectedProperty


def export_sim_db_to_xlsx(jsonFilePath):
    dbFile = open(jsonFilePath)
    dbJson = json.load(dbFile)
    dbDf = pd.json_normalize(dbJson)
    xlsxFilePath = os.path.splitext(jsonFilePath)[0] + '.xlsx'
    dbDf.to_excel(xlsxFilePath)


def generate_fwf(df, formatType='standardDf', outFilePath=None):
    fileTypeSpecs = FILE_TYPES_SPECS[formatType]
    try:
        dfStr = df[fileTypeSpecs['columnOrder']].to_string(
            formatters=fileTypeSpecs['formatters'], index=False,
            header=fileTypeSpecs['header']
        )
    except KeyError as err:
        if err.args[0] in ['columnOrder', 'formatters', 'header']:
            print("FILE_TYPE_SPECS['{:s}']['{:s}'] apparently missing, please verify.".format(
                formatType, err.args[0]
            ))
            raise
        else:
            warnings.warn('No column found.')
            return
    # headerList = standardDf.columns + ' [' + UNITS_STANDARD_DF_EXTENDED + ']'
    if outFilePath is not None:
        outFilePath += FILE_TYPES_SPECS[formatType]['ext']
        with open(outFilePath, 'w') as outFile:
            outFile.write(dfStr)
    return dfStr, outFilePath


def load_standard_fwf(sourceFilePath, removeNanInf=False):
    standardDf = pd.read_fwf(sourceFilePath)
    standardDf = extend_standard_df(standardDf, removeNanInf)
    return standardDf


def convert_irina_distr_to_standard_df(sourceFilePath, filterSpecsSelector=None, outFilePath=None):
    standardDf = pd.read_csv(
        sourceFilePath, delim_whitespace=True, index_col=False,
        header=0, names=('x', 'px', 'y', 'py', 'pz', 't')
    )
    standardDf['z'] = 17.5   # [mm]
    # TODO: Chekc that the units in the following line are correct
    standardDf['t'] = standardDf['t'] / C * 1.e9                                # [ns]
    standardDf['pdgId'] = -11
    standardDf['Q'] = pdgId_to_particle_const(standardDf['pdgId'], 'Q')
    standardDf = extend_standard_df(standardDf)
    if filterSpecsSelector is not None:
        standardDf = use_filter_specs_selector(standardDf, sourceFilePath, filterSpecsSelector)
    if outFilePath is not None:
        _, outFilePath = generate_fwf(standardDf, outFilePath=outFilePath)
    return standardDf, outFilePath


def convert_fcceett_to_standard_df(sourceFilePath, pdgId=[], outFwfPath=None):
    if not isinstance(pdgId, list):
        pdgId = [pdgId]
    rootFile = ROOT.TFile.Open(sourceFilePath, 'READ')
    dfDict = {}
    outFwfPathDict = {}
    for rootKey in rootFile.GetListOfKeys():
        distrName = rootKey.GetName()
        rootTree = rootFile.Get(distrName)
        rootTreeRdf = ROOT.RDataFrame(rootTree)
        rootTreeDict = rootTreeRdf.AsNumpy()
        rootTreeMat = np.column_stack([q for q in rootTreeDict.values()])
        standardDf = pd.DataFrame(data=rootTreeMat, columns=rootTreeDict.keys())
        standardDf.drop('evtId', axis=1, inplace=True)
        standardDf['t'] = standardDf['t'] * 1e-3                                # [ns]
        standardDf['pdgId'] = standardDf['pdgId'].astype(int)
        standardDf['Q'] = pdgId_to_particle_const(standardDf['pdgId'], 'Q')
        # standardDf = extend_standard_df(standardDf)
        # TODO: Check if 'e' in FCC-ee Target Tracking is Ekin or Ekin + Erest (simple curiosity)
        # One could have also used: standardDf.rename(columns={'e': 'Ekin'}, inplace=True)
        fileSuffix = '_' + distrName
        if distrName == 'amor_leave' and pdgId:
            standardDf = standardDf[standardDf['pdgId'].isin(pdgId)]
            fileSuffix += '_pdgId'
            for id in pdgId:
                fileSuffix += '_' + str(id)
        if outFwfPath is not None:
            filePath, fileExt = os.path.splitext(outFwfPath)
            _, tmpOutFwfPath = generate_fwf(standardDf, outFilePath=filePath+fileSuffix+fileExt)
        dfDict[distrName] = standardDf
        outFwfPathDict[distrName] = tmpOutFwfPath
    return dfDict, outFwfPathDict


def convert_astra_to_standard_df(
        sourceFilePath, discardLostParticles=True, filterSpecsSelector=None,
        zProjection=None, zCut=None, removeNanInf=False,
        outFwfPath=None, verbose=False
):
    standardDf = pd.read_csv(
        sourceFilePath, delim_whitespace=True,
        index_col=False, header=None,
        names=FILE_TYPES_SPECS['astra']['columnOrder']
    )
    # names=('x', 'y', 'z', 'px', 'py', 'pz', 't', 'Q', 'pdgId', 'statusFlag')
    standardDf['x'] = standardDf['x'] * 1.e3                                    # [mm]
    standardDf['px'] = standardDf['px'] * 1.e-6                                 # [MeV]
    standardDf['y'] = standardDf['y'] * 1.e3                                    # [mm]
    standardDf['py'] = standardDf['py'] * 1.e-6                                 # [MeV]
    standardDf['z'] = standardDf['z'] * 1.e3                                    # [mm]
    standardDf['pz'] = standardDf['pz'] * 1.e-6                                 # [MeV]
    # Longitudinal coordinates in Astra are with respect to reference particle
    standardDf.loc[1:, 'z'] = standardDf.loc[1:, 'z'] + standardDf.loc[0, 'z']
    standardDf.loc[1:, 'pz'] = standardDf.loc[1:, 'pz'] + standardDf.loc[0, 'pz']
    standardDf.loc[1:, 'clock'] = \
        standardDf.loc[1:, 'clock'] + standardDf.loc[0, 'clock']                  # [ns]
    standardDf.rename(columns={'clock': 't'}, inplace=True)
    standardDf['particleIndex'].replace(to_replace={1: 11, 2: -11}, inplace=True)
    standardDf.rename(columns={'particleIndex': 'pdgId'}, inplace=True)
    standardDf['macroCharge'] = standardDf['macroCharge'] * 1.e-9               # [C]
    standardDf.rename(columns={'macroCharge': 'Q'}, inplace=True)
    if discardLostParticles:
        standardDf = standardDf[standardDf['statusFlag'] >= -6]
    standardDf.drop('statusFlag', axis=1, inplace=True)
    standardDf = extend_standard_df(standardDf, removeNanInf)
    if filterSpecsSelector is not None:
        standardDf = use_filter_specs_selector(standardDf, sourceFilePath, filterSpecsSelector)
    labelStr = ''
    # If requested, cut particles that did not reach zCut
    # TODO: Use function filter_distr?
    if zCut is not None:
        cutInds = standardDf['z'] < zCut
        rejectedParticles = standardDf[cutInds]
        standardDf = standardDf[~cutInds]
        warnings.warn(
            '{:d} particles with z < {:.3f} mm have been discarded.'.format(cutInds.sum(), zCut)
        )
        if verbose and cutInds.sum() > 0:
            print(rejectedParticles)
    # If requested, project all particles on a plane at zProjection
    if zProjection is not None:
        deltaZ = standardDf['z'] - zProjection
        deltaX = deltaZ * standardDf['px'] / standardDf['pz']
        deltaY = deltaZ * standardDf['py'] / standardDf['pz']
        standardDf['x'] = standardDf['x'] - deltaX
        standardDf['y'] = standardDf['y'] - deltaY
        standardDf['z'] = zProjection
        alpha = np.arctan(
            np.sqrt(standardDf['px']**2. + standardDf['py']**2.)
            / standardDf['pz']
        )
        vz = standardDf['betaRel'] * C * np.cos(alpha)
        standardDf['t'] = standardDf['t'] - deltaZ/vz*1e6
        labelStr += '_zProj{:.0f}mm'.format(zProjection)
        # After projection remove inf values generated by pz = 0
        removeNanInf = True
    standardDf = check_nan_inf_in_distr(standardDf, removeNanInf)
    if outFwfPath is not None:
        _, outFwfPath = generate_fwf(standardDf, outFilePath=outFwfPath+labelStr)
    return standardDf, outFwfPath


def convert_rftrack_to_standard_df(
        rftrackDfFormat=None, rftrackDf=None,
        sourceFilePath=None, sourceFormat='octave', octaveMatrixName=None,
        filterSpecsSelector=None, s=np.nan, t=np.nan, pdgId=np.nan, Qbunch=np.nan,
        removeNanInf=False, outFwfPath=None
):
    # TODO: Group input parameters in dictionaries?
    # E.g. {sourceFilePath, sourceFormat, octaveMatrixName}
    if np.isnan(Qbunch):
        warnings.warn(
            'Qbunch has not been declared. ' +
            'The charge Q of the (macro)particles cannot be determined.'
        )
    colNames = FILE_TYPES_SPECS[rftrackDfFormat]['columnOrder']
    if isinstance(rftrackDf, np.ndarray):
        rftrackDf = pd.DataFrame(
            data=rftrackDf,
            columns=colNames
        )
    # TODO: Review use of convert_from_input_check to remove some parameters?
    # This is the only use in a non-convert_standard_df_to_...-function.
    rftrackDf = convert_from_input_check(
        rftrackDf, sourceFilePath, sourceFormat=sourceFormat,
        octaveMatrixName=octaveMatrixName, colNames=colNames
    )
    standardDf = rftrackDf[colNames].copy()
    if rftrackDfFormat == 'rftrack_Px_S':
        # TODO: Just renaming s to z is consistent with the definition of z?
        standardDf.rename(columns={'s': 'z'}, inplace=True)
        standardDf['t'] = t                                                     # [ns]
    elif rftrackDfFormat == 'rftrack_xp_t':
        standardDf['z'] = s                                                     # [mm]
        standardDf['t'] = standardDf['t'] / C * 1e6                             # [ns]
        standardDf = p_components_from_angles(standardDf)
    standardDf['pdgId'] = pdgId
    NparticlesPerBunch = standardDf.shape[0]
    standardDf['Q'] = Qbunch / NparticlesPerBunch                               # [C]
    standardDf = extend_standard_df(standardDf, removeNanInf)
    if filterSpecsSelector is not None:
        standardDf = use_filter_specs_selector(standardDf, sourceFilePath, filterSpecsSelector)
    if outFwfPath is not None:
        _, outFwfPath = generate_fwf(standardDf, outFilePath=outFwfPath)
    return standardDf, outFwfPath


def convert_standard_df_to_rftrack(
        standardDf=None, sourceFilePath=None, rftrackDfFormat=None, outFilePath=None
):
    standardDf = convert_from_input_check(
        standardDf, sourceFilePath
    )
    if rftrackDfFormat == 'rftrack_Px_S':
        rftrackDf = standardDf[['x', 'px', 'y', 'py', 'z']].copy()
        # TODO: Check consistency of definitions, z vs. s
        rftrackDf.rename(columns={'z': 's'}, inplace=True)
    elif rftrackDfFormat == 'rftrack_xp_t':
        rftrackDf = standardDf[['x', 'xp', 'y', 'yp', 't']].copy()
        rftrackDf['t'] = rftrackDf['t'] * C / 1e6  # [mm/c]
    # TODO: Check consistency of definitions, pz vs p
    rftrackDf['p'] = pVect_to_p(standardDf['px'], standardDf['py'], standardDf['pz'])   # [MeV/c]
    if outFilePath is not None:
        _, outFilePath = generate_fwf(
            rftrackDf, formatType=rftrackDfFormat, outFilePath=outFilePath
        )
    return rftrackDf, outFilePath


def convert_standard_df_to_astra(
        standardDf=None, sourceFilePath=None, refParticleId=0,
        outFilePath=None
):
    standardDf = convert_from_input_check(
        standardDf, sourceFilePath
    )
    astraDf = standardDf[['x', 'y', 'z', 'px', 'py', 'pz', 't', 'Q', 'pdgId']].copy()
    astraDf['x'] = astraDf['x'] * 1e-3                                          # [m]
    astraDf['y'] = astraDf['y'] * 1e-3                                          # [m]
    astraDf['z'] = astraDf['z'] * 1e-3                                          # [m]
    astraDf['px'] = astraDf['px'] * 1e6                                         # [eV/c]
    astraDf['py'] = astraDf['py'] * 1e6                                         # [eV/c]
    astraDf['pz'] = astraDf['pz'] * 1e6                                         # [eV/c]
    # t is already in [ns]
    astraDf['Q'] = astraDf['Q'] * 1e9                                           # [nC]
    astraDf['pdgId'].replace(to_replace={11: 1, -11: 2}, inplace=True)
    # Small check for the reference particle
    if not np.isclose(astraDf['x'][refParticleId], 0) \
            or not np.isclose(astraDf['y'][refParticleId], 0):
        warnings.warn(
            'Reference particle is out of axis with (x,y) = ({:.6f}, {:.6f}) m'.format(
                astraDf['x'][refParticleId], astraDf['y'][refParticleId]
            ))
    # Longitudinal variables with respect to the reference particle
    nonRefIds = [ind for ind in standardDf.index if ind != refParticleId]
    astraDf.loc[nonRefIds, 'z'] = astraDf.loc[nonRefIds, 'z'] - astraDf.loc[refParticleId, 'z']
    astraDf.loc[nonRefIds, 'pz'] = astraDf.loc[nonRefIds, 'pz'] - astraDf.loc[refParticleId, 'pz']
    astraDf.loc[nonRefIds, 't'] = astraDf.loc[nonRefIds, 't'] - astraDf.loc[refParticleId, 't']
    # Add status flag
    statusFlag = np.full(astraDf.shape[0], 5)
    astraDf['statusFlag'] = statusFlag
    astraDf['pdgId'] = astraDf['pdgId'].astype(int)
    astraDf['statusFlag'] = astraDf['statusFlag'].astype(int)
    astraDf = astraDf.loc[[refParticleId]+nonRefIds, :]
    astraDf.columns = FILE_TYPES_SPECS['astra']['columnOrder']
    if outFilePath is not None:
        _, outFilePath = generate_fwf(astraDf, formatType='astra', outFilePath=outFilePath)
    return astraDf, outFilePath


def convert_placet_to_standard_df(
        sourceFilePath, filterSpecsSelector=None, z0=0, pdgId=-11, Qbunch=np.nan,
        removeNanInf=False, outFwfPath=None
):
    if np.isnan(Qbunch):
        warnings.warn(
            'Qbunch has not been declared. ' +
            'The charge Q of the (macro)particles cannot be determined.'
        )
    standardDf = pd.read_csv(
        sourceFilePath, engine='python', delim_whitespace=True,
        index_col=False, header=None,
        names=FILE_TYPES_SPECS['placet']['columnOrder'],
        skiprows=5
    )
    standardDf['x'] = standardDf['x'] * 1e-3  # From [um] to [mm]
    standardDf['y'] = standardDf['y'] * 1e-3  # From [um] to [mm]
    standardDf['p'] = standardDf['p'] * 1e3   # From [GeV/c] to [MeV/c]
    standardDf['xp'] = standardDf['xp'] * 1e-3  # From [urad] to [mrad]
    standardDf['yp'] = standardDf['yp'] * 1e-3  # From [urad] to [mrad]
    standardDf = p_components_from_angles(standardDf)
    standardDf['z'] = z0  # [mm]
    standardDf['t'] = standardDf['t'] / C * 1e3  # From [um/c] to [ns]
    standardDf['pdgId'] = pdgId
    NparticlesPerBunch = standardDf.shape[0]
    standardDf['Q'] = Qbunch / NparticlesPerBunch  # [C]
    standardDf = extend_standard_df(standardDf, removeNanInf)
    if filterSpecsSelector is not None:
        standardDf = use_filter_specs_selector(standardDf, sourceFilePath, filterSpecsSelector)
    if outFwfPath is not None:
        _, outFwfPath = generate_fwf(standardDf, outFilePath=outFwfPath)
    return standardDf, outFwfPath


def convert_standard_df_to_placet(
        standardDf=None, sourceFilePath=None, outFilePath=None
):
    standardDf = convert_from_input_check(
        standardDf, sourceFilePath
    )
    placetDf = standardDf[['x', 'y', 't', 'xp', 'yp']].copy()
    placetDf['x'] = placetDf['x'] * 1e3  # [um]
    placetDf['y'] = placetDf['y'] * 1e3  # [um]
    placetDf['t'] = placetDf['t'] * C * 1e3  # [um/c]
    placetDf['xp'] = placetDf['xp'] * 1e3  # [urad]
    placetDf['yp'] = placetDf['yp'] * 1e3  # [urad]
    placetDf['p'] = pVect_to_p(
        standardDf['px'], standardDf['py'], standardDf['pz']
    ) * 1e-3  # [GeV/c]
    if outFilePath is not None:
        _, outFilePath = generate_fwf(placetDf, formatType='placet', outFilePath=outFilePath)
    return placetDf, outFilePath


def convert_from_input_check(
        df, sourceFilePath, sourceFormat='standardDf',
        octaveMatrixName=None, colNames=None
):
    if sourceFilePath is not None:
        if df is None:
            if sourceFormat == 'standardDf':
                df = load_standard_fwf(sourceFilePath)
            elif sourceFormat == 'octave':
                df = opi.load_octave_matrices(
                    sourceFilePath, matNamesToLoad=octaveMatrixName,
                    colNames=colNames
                )
        else:
            raise ValueError(
                'Only one between df and sourceFilePath can be set different than None.'
            )
    elif df is None:
        raise ValueError(
            'At least one between df and sourceFilePath must be different than None.'
        )
    return df


def convert_sdds_to_standard_df(
        sourceFilePath, filterSpecsSelector=None, z0=0, pdgId=-11, Qbunch=np.nan,
        removeNanInf=False, outFwfPath=None
):
    os.system('sddsconvert -ascii ' + sourceFilePath)
    standardDf = pd.read_csv(
        sourceFilePath, skiprows=24, delim_whitespace=True,
        names=FILE_TYPES_SPECS['sdds']['columnOrder']
    )
    standardDf.drop('trackingId', axis=1, inplace=True)
    pCentral = get_sdds_parameter(sourceFilePath, 'pCentral')
    QbunchRead = get_sdds_parameter(sourceFilePath, 'Charge')
    if not np.isnan(Qbunch) and not np.isclose(Qbunch, QbunchRead):
        warnStr = 'Declared Qbunch = {:.6g} C (function input) '.format(Qbunch) + \
            'does not correspond to value in file QbunchRead = {:.6f} C. '.format(QbunchRead) + \
            'Going to use '
        if not np.isnan(QbunchRead):
            Qbunch = QbunchRead
            warnStr += 'QbunchRead.'
        else:
            warnStr += 'declared Qbunch (function input).'
        warnings.warn(warnStr)
    else:
        Qbunch = QbunchRead
    NparticlesPerBunch = get_sdds_parameter(sourceFilePath, 'Particles')
    Erest = pdgId_to_particle_const(standardDf['pdgId'], 'Erest')
    standardDf['pdgId'] = pdgId
    standardDf['x'] = standardDf['x'] * 1.e3                                    # [mm]
    standardDf['y'] = standardDf['y'] * 1.e3                                    # [mm]
    standardDf['z'] = z0                                                        # [mm]
    standardDf['t'] = standardDf['t'] * 1.e9                                    # [ns]
    standardDf['xp'] = standardDf['xp'] * 1.e3                                  # [mrad]
    standardDf['yp'] = standardDf['yp'] * 1.e3                                  # [mrad]
    standardDf['p'] = standardDf['p'] * Erest                                   # [MeV/c]
    standardDf = p_components_from_angles(standardDf)
    standardDf['Q'] = Qbunch / NparticlesPerBunch                               # [C]
    # Some extended variables already available
    standardDf = extend_standard_df(standardDf, removeNanInf)
    # sElegant = standardDf['betaRel'] * C * standardDf['t'] * 1e-6             # [mm]
    if filterSpecsSelector is not None:
        standardDf = use_filter_specs_selector(standardDf, sourceFilePath, filterSpecsSelector)
    if outFwfPath is not None:
        _, outFwfPath = generate_fwf(standardDf, outFilePath=outFwfPath)
    return standardDf, outFwfPath


def p_components_from_angles(standardDf):
    p = standardDf['p']  # [MeV/c]
    standardDf.drop('p', axis=1, inplace=True)
    standardDf['pz'] = p / np.sqrt(
        1. + np.tan(standardDf['xp']*1e-3)**2. + np.tan(standardDf['yp']*1e-3)**2.
    )  # [MeV/c]
    halfPi = np.pi / 2. * 1e3  # [mrad]
    backwardsMovingInds = (
            ((standardDf['xp'] > halfPi) | (standardDf['xp'] < -halfPi))
            & ((standardDf['yp'] >= -halfPi) & (standardDf['yp'] <= halfPi))
        ) | (
            ((standardDf['yp'] > halfPi) | (standardDf['yp'] < -halfPi))
            & ((standardDf['xp'] >= -halfPi) & (standardDf['xp'] <= halfPi))
        )
    # TODO: Check what happens when both angles are larger than pi/2.
    standardDf.loc[backwardsMovingInds, 'pz'] *= -1
    standardDf['px'] = np.tan(standardDf['xp']*1e-3) * standardDf['pz']  # [MeV/c]
    standardDf['py'] = np.tan(standardDf['yp']*1e-3) * standardDf['pz']  # [MeV/c]
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
        warnings.warn("Could not read parameter {:s} correctly (resultStr = '{:s}').".format(
            parameterName, resultStr
        ))
    return resultVal


def convert_standard_df_to_sdds(
        standardDf=None, sourceFilePath=None, refParticleId=0, outFilePath=None
):
    standardDf = convert_from_input_check(
        standardDf, sourceFilePath
    )
    if outFilePath is None:
        raise ValueError('Please specify outFilePath.')
    astraDf, _ = convert_standard_df_to_astra(
        standardDf=standardDf, refParticleId=refParticleId
    )
    astraDfStr, _ = generate_fwf(astraDf, formatType='astra', outFilePath='TmpAstra')
    outFilePath += FILE_TYPES_SPECS['sdds']['ext']
    os.system('astra2elegant TmpAstra.001 ' + outFilePath)
    # os.system('rm TmpAstra.001')
    # os.system('astra2elegant -pipe=in ' + astraDfStr + ' ' + outFilePath)


def check_nan_inf_in_distr(distrIn, removeNanInf):
    NpartIn = distrIn.shape[0]
    # TODO: Better define subset
    colsToCheck = ['x', 'y', 'px', 'py', 'xp', 'yp']
    distrWithoutNan = distrIn.dropna(subset=colsToCheck)
    NpartWithoutNan = distrWithoutNan.shape[0]
    if removeNanInf:
        strMessage = 'removed'
    else:
        strMessage = 'found'
    if NpartIn - NpartWithoutNan > 0:
        warnings.warn('{:d} particles with some NaN value have been {:s}:'.format(
            NpartIn - NpartWithoutNan, strMessage
        ))
        print(distrIn[distrIn[colsToCheck].isna().any(axis=1)])
    # TODO: Better define subset
    distrWithoutNanInf = distrWithoutNan.replace([np.inf, -np.inf], np.nan) \
        .dropna(subset=colsToCheck)
    NpartWithoutNanInf = distrWithoutNanInf.shape[0]
    if NpartWithoutNan - NpartWithoutNanInf > 0:
        warnings.warn('{:d} particles with some Inf value have been {:s}:'.format(
            NpartWithoutNan - NpartWithoutNanInf, strMessage
        ))
        print(distrIn[distrIn[colsToCheck].isin([np.inf, -np.inf]).any(axis=1)])
    if removeNanInf:
        return distrWithoutNanInf
    else:
        return distrIn


def plot_hist(
        ax, distr, binWidth=None, binLims=None, density=False, legendLabel='',
        orientation='vertical', parsInLabel=True, opacityHist=1.
):
    defaultBinNum = 100
    if binLims is None:
        if binWidth is None:
            binWidth = (np.nanmax(distr) - np.nanmin(distr)) / (defaultBinNum-1)
            if binWidth == 0:
                binWidth = 1.
                defaultBinNum = 3
        binLims = [np.nanmin(distr)-binWidth*1.5, np.nanmax(distr)+binWidth*1.5]
    if binWidth is None:
        binEdges = np.linspace(binLims[0], binLims[1], num=defaultBinNum)
    else:
        binEdges = np.arange(binLims[0], binLims[1], binWidth)
    counts, _, histObj = ax.hist(
        distr, bins=binEdges, density=density,
        orientation=orientation, alpha=opacityHist
    )
    # Select portion of distribution
    distr = distr[(distr >= np.nanmin(binEdges)) & (distr <= np.nanmax(binEdges))]
    avg = np.mean(distr)
    std = np.std(distr)
    # TODO: Allows for input of precision of label values
    if parsInLabel:
        if legendLabel != '':
            legendLabel += ', '
        legendLabel += 'avg = {:.2e}, std = {:.2e}'.format(avg, std)
    # Plot Gaussian fit
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        try:
            mu, sigma = scistats.norm.fit(distr)
            binCenters = (binEdges[1:] + binEdges[:-1]) / 2.
            gaussFits = scistats.norm.pdf(binCenters, mu, sigma)
        except RuntimeWarning as err:
            knownMessages = [
                'divide by zero encountered in divide',
                'divide by zero encountered in true_divide',
                'invalid value encountered in divide',
                'invalid value encountered in true_divide'
            ]
            if not str(err) in knownMessages:
                raise
        else:
            binWidths = binEdges[1:] - binEdges[:-1]
            fittedHist = gaussFits * binWidths * distr.shape[0]
            if orientation == 'vertical':
                data = [binCenters, fittedHist]
            elif orientation == 'horizontal':
                data = [fittedHist, binCenters]
            ax.plot(
                *data, '-', color=histObj.get_children()[0].get_facecolor()
            )
            # TODO: Why is sigmaGauss sometimes different than std?
            # if parsInLabel:
            #     legendLabel += ', sigmaGauss = {:.2e}'.format(sigma)
    histObj.set_label(legendLabel)
    ax.legend()
    return avg, std


def set_lims(ax, direction, var, lims):
    if lims is None:
        lims = (np.nanmin(var), np.nanmax(var))
    if not np.isclose(*lims):
        if direction == 'x':
            ax.set_xlim(lims)
        elif direction == 'y':
            ax.set_ylim(lims)


def scatter_individual_marker_style(pathCollectionObj, markerStyles):
    if markerStyles is None:
        return
    elif isinstance(markerStyles, str):
        markerStyles = [markerStyles] * pathCollectionObj.get_transforms().shape[0]
    markerPaths = []
    for m in markerStyles:
        markerObj = pltMarkers.MarkerStyle(m, 'none')
        markerPath = markerObj.get_path().transformed(markerObj.get_transform())
        markerPaths.append(markerPath)
    pathCollectionObj.set_paths(markerPaths)
    # Next line: Explicitly set the edge color before removing the facecolor
    pathCollectionObj.set_edgecolor(pathCollectionObj.get_facecolor())
    pathCollectionObj.set_facecolor('none')


def plot_phase_space_2d(
        ax, distr, varName1=None, varName2=None, var1Center=False, var2Center=False,
        title=None, legendLabel='', displayTable=True,
        binWidth1=None, binWidth2=None, lims1=None, lims2=None, pzCutoff=None,
        markerStyle=None, markerSize=2, color=None, opacityHist=1.
):
    if varName1 is None or varName2 is None:
        raise ValueError('varName1 and varName2 are mandatory arguments.')
    opacityScatter = 1.
    x = None
    y = None
    for ind, (vName, vCenter, var) in enumerate(
            [[varName1, var1Center, x], [varName2, var2Center, y]]):
        try:
            var = distr[vName]
        except KeyError:
            if vName == 'r':
                var = np.sqrt(distr['x']**2. + distr['y']**2.)
        if vCenter:
            # TODO: Does this change of distr also remains outside of the function?
            var = var - var.mean()
        if ind == 0:
            x = var
        else:
            y = var
    if color is None:
        # Calculate the point density
        # xy = np.vstack([distr[varName1], distr[varName2]])
        # color = scistats.gaussian_kde(xy)(xy)
        # TODO: Check number of bins.
        # Are small density region not plotted if bins is small (e.g. 100)?
        data, x_e, y_e = np.histogram2d(x, y, bins=1000, density=True)
        try:
            color = sciinterp.interpn(
                (0.5*(x_e[1:]+x_e[:-1]), 0.5*(y_e[1:]+y_e[:-1])), data, np.vstack([x, y]).T,
                method='splinef2d', bounds_error=False
            )
            if np.all(np.isnan(color)):
                color = np.ones(x.shape)
        except ValueError:
            color = np.ones(x.shape)
        # Sort the points by density, so that the densest points are plotted last
        idx = color.argsort()
        x, y, color = x.iloc[idx], y.iloc[idx], color[idx]
    pathColl = ax[0, 0].scatter(
        x, y,
        s=markerSize, c=color, marker='.', alpha=opacityScatter
    )
    scatter_individual_marker_style(pathColl, markerStyle)
    plot_hist(
        ax[1, 0], x, binWidth=binWidth1, binLims=lims1,
        legendLabel=legendLabel, opacityHist=opacityHist
    )
    plot_hist(
        ax[0, 1], y, binWidth=binWidth2, binLims=lims2,
        legendLabel=legendLabel, orientation='horizontal',
        opacityHist=opacityHist
    )
    pzLabels = ['All, Counts = {:d}'.format(x.shape[0]), ]
    if pzCutoff is not None:
        distrLowPz = distr[distr['pz'] < pzCutoff]
        ax[0, 0].scatter(
            distrLowPz[varName1], distrLowPz[varName2],
            s=markerSize, marker=markerStyle, facecolors='none', linewidths=1., edgecolors=color,
            alpha=opacityScatter
        )
        plot_hist(
            ax[1, 0], distrLowPz[varName1],
            binWidth=binWidth1, binLims=lims1, opacityHist=opacityHist
        )
        plot_hist(
            ax[0, 1], distrLowPz[varName2], orientation='horizontal',
            binWidth=binWidth2, binLims=lims2, opacityHist=opacityHist
        )
        pzLabels += ['pz <= {:.1f} {:s}, Counts = {:d}'.format(
            pzCutoff, UNITS_STANDARD_DF_EXTENDED['pz'], distrLowPz.shape[0]
        )]
        ax[0, 0].legend(
            pzLabels, markerscale=5.,
            loc='upper left', bbox_to_anchor=(1.2, -0.2)
        )
    set_lims(ax[0, 0], 'x', x, lims1)
    set_lims(ax[0, 0], 'y', y, lims2)
    ax[1, 0].set_xlim(ax[0, 0].get_xlim())
    ax[1, 0].invert_yaxis()
    # ax[1,0].set_xticklabels([])
    ax[0, 1].set_ylim(ax[0, 0].get_ylim())
    # ax[0,1].set_yticklabels([])
    ax[0, 0].grid(True)
    ax[0, 1].grid(True)
    ax[1, 0].grid(True)
    ax[0, 0].set_xlabel(varName1+' ['+UNITS_STANDARD_DF_EXTENDED[varName1]+']')
    ax[0, 0].set_ylabel(varName2+' ['+UNITS_STANDARD_DF_EXTENDED[varName2]+']')
    # ax[0,0].yaxis.set_label_position('right')
    # ax[0,0].yaxis.tick_right()
    ax[1, 0].set_ylabel('Bin counts')
    ax[0, 1].set_xlabel('Bin counts')
    if displayTable and (varName1 == 'x' and varName2 == 'xp') \
            or (varName1 == 'y' and varName2 == 'yp'):
        # TODO: Pass filter specs
        plot_parameters(ax[1, 1], distr, varName1, filterSpecs={})
    else:
        ax[1, 1].set_visible(False)
    if title is not None:
        ax[0, 0].get_figure().suptitle(title)


def plot_parameters(ax, distr, planeName, filterSpecs={}):
    emitTab = []
    emitDefinitions = ['normalized', 'geometric', 'tracespace']
    # vPos = 0.1
    # alignSpecs = {
    #     'horizontalalignment': 'left',
    #     'verticalalignment': 'center',
    #     'transform': ax.transAxes
    # }
    for emitDef in emitDefinitions:
        emitTab.append(compute_emittance(
            distr, planeName, norm=emitDef,
            filterSpecs=filterSpecs, verbose=False
        ))
    emitDefinitions = [
        'emit'+emitDef.capitalize() for emitDef in emitDefinitions
    ]
    # TODO: Parametrize emittance units
    # TODO: EMIT_DEFINITIONS and TWISS_NAMES
    twissNames = ['alphaTwiss', 'betaTwiss', 'gammaTwiss']
    emitTab += compute_twiss(
        distr, planeName, filterSpecs=filterSpecs, verbose=False
    )
    emitDf = pd.DataFrame([emitTab, ], columns=emitDefinitions+twissNames)
    emitDf.update(emitDf.applymap('{:.3f}'.format))
    tab = pd.plotting.table(
        ax, emitDf.T, loc='center',
        colWidths=[1./(1.+emitDf.shape[0])] * emitDf.shape[0]
    )
    tab = tab.get_celld()
    for cell in tab.values():
        cell.set_height(0.1)
    ax.set_axis_off()


def plot_distr(
        distributions, plotDefs, colors=None, markerStyle=None, markerSize=15,
        title=None, legendLabels=None, figHeight=6.4, figWidth=9.6
):
    # TODO: Integrate handling of distributions and legendLabels in check_marker_specs?
    if isinstance(distributions, pd.DataFrame):
        distributions = [distributions]
        if isinstance(legendLabels, str):
            legendLabels = [legendLabels]
    if legendLabels is not None and len(distributions) != len(legendLabels):
        raise ValueError(
            'len(distributions) = {:d} does not match len(legendLabels) = {:d}.'.format(
                len(distributions), len(legendLabels)
            ))
    if legendLabels is None:
        legendLabels = [''] * len(distributions)
    markerStyle = check_marker_specs(distributions, markerStyle, 'markerStyle')
    markerSize = check_marker_specs(distributions, markerSize, 'markerSize')
    if colors is None:
        if len(distributions) == 1:
            colors = [None]
        else:
            # TODO: Place default colors at the beginning of the module?
            defaultColors = plt.rcParams['axes.prop_cycle'].by_key()['color']
            colors = defaultColors[:len(distributions)]
    axList = []
    for plotDef in plotDefs:
        fig, ax = plt.subplots(2, 2, figsize=(figWidth, figHeight))
        distrAndProperties = zip(
            distributions, markerStyle, markerSize, legendLabels, colors
        )
        for distr, mStyle, mSize, label, color in distrAndProperties:
            if len(distributions) > 1:
                plotDef['pzCutoff'] = None
            plot_phase_space_2d(
                ax, distr, **plotDef,
                markerStyle=mStyle, markerSize=mSize, color=color,
                title=title, legendLabel=label
            )
        axList.append(ax)
    return axList


def plot_ellipse(
        ax, emitGeom, semiAxisOrder=1, color='k',
        alphaTwiss=None, betaTwiss=None, gammaTwiss=None
):
    alphaTwiss, betaTwiss, gammaTwiss = third_twiss_param(
        alphaTwiss=alphaTwiss, betaTwiss=betaTwiss
    )
    tiltAngle = np.rad2deg(
        np.arctan(2*alphaTwiss/(gammaTwiss-betaTwiss)) / 2.
    )
    sqrtTerm = np.sqrt((betaTwiss+gammaTwiss)**2. - 4.)
    majorSemiAxis = 1. / np.sqrt((betaTwiss+gammaTwiss+sqrtTerm)/2./emitGeom)
    minorSemiAxis = 1. / np.sqrt((betaTwiss+gammaTwiss-sqrtTerm)/2./emitGeom)
    if semiAxisOrder == 1:
        semiAxis1 = majorSemiAxis
        semiAxis2 = minorSemiAxis
    elif semiAxisOrder == 2:
        semiAxis1 = minorSemiAxis
        semiAxis2 = majorSemiAxis
    ellipse = pltPatches.Ellipse(
        (0, 0), 2.*semiAxis1, 2.*semiAxis2, tiltAngle,
        facecolor='none', edgecolor=color
    )
    ax.add_patch(ellipse)


def third_twiss_param(alphaTwiss=None, betaTwiss=None, gammaTwiss=None):
    if alphaTwiss is not None and betaTwiss is not None and gammaTwiss is None:
        gammaTwiss = (1. + alphaTwiss**2.) / betaTwiss
    elif alphaTwiss is not None and gammaTwiss is not None and betaTwiss is None:
        betaTwiss = (1. + alphaTwiss**2.) / gammaTwiss
    elif betaTwiss is not None and gammaTwiss is not None and alphaTwiss is None:
        alphaTwiss = np.sqrt(betaTwiss*gammaTwiss - 1.)
    else:
        raise ValueError(
            'Wrong input. Two out of three Twiss parameters are required.'
        )
    return alphaTwiss, betaTwiss, gammaTwiss


def distr_within_ellipse(standardDf, emitTraceSpace, ellipseSpecs):
    indsWithinEllipse = pd.Series(True, index=standardDf.index)
    for planeName, ellSpecs in ellipseSpecs.items():
        alphaTwiss, betaTwiss, gammaTwiss = third_twiss_param(**ellSpecs)
        u = standardDf[planeName]
        uDiv = standardDf[planeName+'p']
        # TODO: Allow for different emittances in the different planes
        indsWithinEllipse = indsWithinEllipse & (
            gammaTwiss*u**2. + 2.*alphaTwiss*u*uDiv + betaTwiss*uDiv**2.
            <= emitTraceSpace
        )
    distrWithinEllipse = standardDf[indsWithinEllipse]
    portion = distrWithinEllipse.shape[0] / standardDf.shape[0]
    return distrWithinEllipse, portion


def check_marker_specs(distributions, markerSpecs, specName):
    try:
        if len(markerSpecs) != len(distributions):
            raise ValueError('{:s} and distributions must have the same length.'.format(specName))
    except TypeError:
        markerSpecs = [markerSpecs] * len(distributions)
    return markerSpecs


def set_plot_defs_from_distrs(distrList, setNames):
    distrAll = pd.concat(distrList)
    distrAll = check_nan_inf_in_distr(distrAll, True)
    plotDefs = {}
    plotDefs['TransvPlane'] = [
        {
            'varName1': 'x', 'varName2': 'y',
            'lims1': (distrAll['x'].min(), distrAll['x'].max()),
            'lims2': (distrAll['y'].min(), distrAll['y'].max()),
            'opacityHist': 0.6,
        }
    ]
    plotDefs['TransvPsAngles'] = [
        {
            'varName1': 'x', 'varName2': 'xp',
            'lims1': (distrAll['x'].min(), distrAll['x'].max()),
            'lims2': (distrAll['xp'].min(), distrAll['xp'].max()),
            'opacityHist': 0.6,
        },
        {
            'varName1': 'y', 'varName2': 'yp',
            'lims1': (distrAll['y'].min(), distrAll['y'].max()),
            'lims2': (distrAll['yp'].min(), distrAll['yp'].max()),
            'opacityHist': 0.6,
        }
    ]
    plotDefs['LongPsZ'] = [
        {
            'varName1': 'z', 'varName2': 'pz',
            'lims1': (distrAll['z'].min(), distrAll['z'].max()),
            'lims2': (distrAll['pz'].min(), distrAll['pz'].max()),
            'opacityHist': 0.6,
        },
        {
            'varName1': 'z', 'varName2': 'Ekin',
            'lims1': (distrAll['z'].min(), distrAll['z'].max()),
            'lims2': (distrAll['Ekin'].min(), distrAll['Ekin'].max()),
            'opacityHist': 0.6,
        }
    ]
    plotDefs['LongPsT'] = [
        {
            'varName1': 't', 'varName2': 'pz',
            'lims1': (distrAll['t'].min(), distrAll['t'].max()),
            'lims2': (distrAll['pz'].min(), distrAll['pz'].max()),
            'opacityHist': 0.6,
        },
        {
            'varName1': 't', 'varName2': 'Ekin',
            'lims1': (distrAll['t'].min(), distrAll['t'].max()),
            'lims2': (distrAll['Ekin'].min(), distrAll['Ekin'].max()),
            'opacityHist': 0.6,
        }
    ]
    plotDefs['CoordsVsT'] = [
        {
            'varName1': 't', 'varName2': 'x',
            'lims1': (distrAll['t'].min(), distrAll['t'].max()),
            'lims2': (distrAll['x'].min(), distrAll['x'].max()),
            'opacityHist': 0.6,
        },
        {
            'varName1': 't', 'varName2': 'y',
            'lims1': (distrAll['t'].min(), distrAll['t'].max()),
            'lims2': (distrAll['y'].min(), distrAll['y'].max()),
            'opacityHist': 0.6,
        },
        {
            'varName1': 't', 'varName2': 'z',
            'lims1': (distrAll['t'].min(), distrAll['t'].max()),
            'lims2': (distrAll['z'].min(), distrAll['z'].max()),
            'opacityHist': 0.6,
        },
    ]
    plotDefs['r-pz'] = [
        {
            'varName1': 'r', 'varName2': 'pz',
            'lims1': (0, np.sqrt((distrAll['x']**2. + distrAll['y']**2.).max())),
            'lims2': (distrAll['pz'].min(), distrAll['pz'].max()),
            'opacityHist': 0.6,
        },
    ]
    plotDefsCollection = []
    for setName in setNames:
        plotDefsCollection += plotDefs[setName]
    return plotDefsCollection


def generate_fieldmap_astra_ideal_tw(fileBasePath, freq, Lstructure, zRes):
    z = np.arange(0, Lstructure+zRes, zRes)
    # k = 2*np.pi*freq / C
    # Ereal = np.cos(k*z)
    # Eimag = np.sin(k*z)
    # np.savetxt(fileBasePath+'_Real.astra', np.stack([z, Ereal], axis=1), fmt='%.9e')
    # np.savetxt(fileBasePath+'_Imag.astra', np.stack([z, Eimag], axis=1), fmt='%.9e')
    Eabs = np.array([1.0]*len(z))
    np.savetxt(fileBasePath+'.astra', np.stack([z, Eabs], axis=1), fmt='%.9e')


def generate_lattice_quad_over_rf_elegant(
        elementName, L, Nslices, quadGradient, rfFreq, rfPhase, rfVoltage, EkinIni,
        pdgId=-11, ZoverA=1., elegantInputFilePath=None, quadOrder=2
):
    inputStr = ''
    inputStr += 'Q: CHARGE, TOTAL=5.0e-9\n'
    Lslice = L / Nslices
    rfVoltageSlice = rfVoltage / Nslices
    # RF
    inputStr += '{:s}.Rf.HalfSlice: RFCA, FREQ={:.1f}, PHASE={:.3f}, &\n'.format(
        elementName, rfFreq, rfPhase
    )
    inputStr += '\tL={:.6f}, VOLT={:.5e}, CHANGE_P0=1, &\n'.format(
        Lslice/2., rfVoltageSlice*1e6/2.
    )
    inputStr += '\tEND1_FOCUS=0, END2_FOCUS=0, BODY_FOCUS_MODEL="NONE"\n'
    # ND
    inputStr += '{:s}.Nd.HalfSlice: DRIFT, L={:.6f}\n'.format(
        elementName, -Lslice/2.
    )
    quadOverRfLineStr = ''
    p0SliceStart = Ekin_to_p(EkinIni, pdgId)
    for sliceInd in range(1, Nslices+1):
        beta = Ekin_to_beta(EkinIni, -11)
        quadStrengthSlice = quad_strength(
            quadGradient, p0SliceStart+rfVoltageSlice/2., ZoverA=1.
        )
        # QUADRUPOLE
        inputStr += '{:s}.Quad.{:d}: QUADRUPOLE, L={:.6f}, K1={:.5e}, ORDER={:d}\n'.format(
            elementName, sliceInd, Lslice, quadStrengthSlice, quadOrder
        )
        inputStr += (
            '{0:s}.Slice{1:d}: Line=({0:s}.Rf.HalfSlice, {0:s}.Nd.HalfSlice, ' +
            '{0:s}.Quad.{1:d}, {0:s}.Nd.HalfSlice, {0:s}.Rf.HalfSlice)\n'
        ).format(elementName, sliceInd)
        # inputStr += generate_quad_over_rf_slice(
        #     Lslice, sliceInd, quadStrengthSlice, rfFreq, rfPhase, rfVoltageSlice
        # )
        quadOverRfLineStr += '{:s}.Slice{:d}'.format(elementName, sliceInd)
        if sliceInd < Nslices:
            quadOverRfLineStr += ', '
        p0SliceStart += rfVoltageSlice
    inputStr += '{:s}.Line: Line=({:s})\n'.format(
        elementName, quadOverRfLineStr
    )
    inputStr += '{:s}.Start: MARKER\n'.format(elementName)
    inputStr += '{:s}.End: MARKER\n'.format(elementName)
    inputStr += '{0:s}: Line=(Q, {0:s}.Start, {0:s}.Line, {0:s}.End)\n'.format(
        elementName
    )
    if elegantInputFilePath is not None:
        with open(elegantInputFilePath, 'w') as outFile:
            outFile.write(inputStr)
    return inputStr


def quad_strength(quadGradient, p, ZoverA=1.):
    quadStrength = C * 1e-6 * ZoverA * quadGradient / p
    return quadStrength


def quad_gradient(quadStrength, p, ZoverA=1.):
    quadGradient = quadStrength * p / C / 1e-6 / ZoverA
    return quadGradient


def quad_matrix(k, Lquad):
    sqrtK = np.sqrt(np.abs(k))
    phi = sqrtK * Lquad
    MquadFocusing = np.array((
        (np.cos(phi), np.sin(phi)/sqrtK),
        (-np.sin(phi)*sqrtK, np.cos(phi))
    ))
    MquadDefocusing = np.array((
        (np.cosh(phi), np.sinh(phi)/sqrtK),
        (np.sinh(phi)*sqrtK, np.cosh(phi))
    ))
    return MquadFocusing, MquadDefocusing


def generate_cross_distribution(
        xMax, yMax, p0, pzDelta, xPoints=5, yPoints=5, pzPoints=5, outFilePath=None
):
    if not isinstance(xPoints, int) or xPoints < 1:
        raise ValueError('xPoints must be an integer >= 1.')
    if not isinstance(yPoints, int) or yPoints < 1:
        raise ValueError('yPoints must be an integer >= 1.')
    if not isinstance(pzPoints, int) or pzPoints < 1:
        raise ValueError('pzPoints must be an integer >= 1.')
    # Define variations
    xArray = np.linspace(-xMax, xMax, num=xPoints)
    yArray = np.linspace(-yMax, yMax, num=yPoints)
    pzArray = p0 + np.linspace(-pzDelta, pzDelta, num=pzPoints)
    # Prepare all combinations of the variations
    xArray = np.repeat(xArray, pzPoints)
    yArray = np.repeat(yArray, pzPoints)
    xArray = np.concatenate((xArray, np.zeros(yPoints*pzPoints)))
    yArray = np.concatenate((np.zeros(xPoints*pzPoints), yArray))
    pzArray = np.tile(pzArray, xPoints+yPoints)
    # Define remaining parameters
    totParticles = (xPoints + yPoints) * pzPoints
    zArray = np.zeros(totParticles)
    pxArray = np.zeros(totParticles)
    pyArray = np.zeros(totParticles)
    tArray = np.zeros(totParticles)
    pdgId = -11
    pdgIdArray = np.full(totParticles, pdgId)
    QArray = pdgId_to_particle_const(pdgIdArray, 'Q')
    dfData = np.stack((
        xArray, pxArray, yArray, pyArray, zArray, pzArray,
        tArray, pdgIdArray, QArray
    ), axis=1)
    standardDf = pd.DataFrame(
        dfData, columns=FILE_TYPES_SPECS['standardDf']['columnOrder']
    )
    standardDf.drop_duplicates(ignore_index=True, inplace=True)
    standardDf = extend_standard_df(standardDf, False)
    if outFilePath is not None:
        _, outFilePath = generate_fwf(standardDf, outFilePath=outFilePath)
    return standardDf, outFilePath


def generate_solenoid_fieldmap(L, B0, Rin):
    """Create solenoid fieldmap.

    Parameters
    ----------
    L : :obj:`float`
        Solenoid length in [m].
    Bz : :obj:`float`
        Solenoidal field in [T].
    Rin : :obj:`float`
        Solenoid internal radius (aperture) in [m].

    Returns
    -------
    zAxis : :obj:`float`
        Sampling positions z of the fieldmap in [m].
    BzOnAxis : :obj:`numpy.array`
        Longitudinal magnetic field on axis in [T].

    """
    zAxis = np.linspace(-2*L, 2*L, 200)
    lp = 0.5 * L + zAxis   # [m]
    lm = 0.5 * L - zAxis   # [m]
    BzOnAxis = 0.5*B0*(lm / np.hypot(Rin, lm) + lp / np.hypot(Rin, lp))   # [T]
    n_L = 95 / 0.236   # [turns/m]  (CLEAR solenoid)
    # T / mu0 / (1/m) = 795774.7150238460162655 A
    print(f'Solenoid current = {B0/n_L*795774.7150238460162655} A\n')
    return zAxis, BzOnAxis


def generate_solenoid_fieldmap_wilson(zAxis, zCenter, RinCoil, RoutCoil, LhalfCoil, J):
    """Create solenoid fieldmap, formula from Wilson, chapter 3, function from Jaap Kosse (PSI).

    Parameters
    ----------
    zAxis : :obj:`float`
        Sampling positions z of the fieldmap in [m].
    zCenter : :obj:`float`
        Position z of the coil center in [m].
    RinCoil : :obj:`float`
        Internal radius of the coil [m].
    RoutCoil : :obj:`float`
        Outer radius of the coil [m].
    LhalfCoil : :obj:`float`
        Half length of the coil [m].
    J : :obj:`float`
        Current density [A/m^2].

    Returns
    -------
    zAxis : :obj:`float`
        Sampling positions z of the fieldmap in [m].
    BzOnAxis : :obj:`numpy.array`
        Longitudinal magnetic field on axis in [T].

    """
    alpha = RoutCoil / RinCoil
    beta1 = (LhalfCoil+zCenter - zAxis) / RinCoil
    beta2 = (LhalfCoil-zCenter + zAxis) / RinCoil
    BzOnAxis = 0.5 * J * RinCoil * (f_factor(alpha, beta1)+f_factor(alpha, beta2))
    return BzOnAxis


def f_factor(alpha, beta):
    """Compute F factor, help function to create solenoid fieldmap.

    Formula from Wilson, chapter 3, function from Jaap Kosse (PSI).

    Parameters
    ----------
    alpha : :obj:`float`
        Parameter alpha (dimensionless).
    beta : :obj:`float`
        Parameter beta (dimensionless).

    Returns
    -------
    F : :obj:`float`
        F factor [N/A^2].

    """
    mu0 = 4 * np.pi * 1e-7   # Vacuum permeability [H/m = N/A^2]
    F = mu0 * beta * np.log((alpha+(alpha**2.+beta**2.)**0.5) / (1.+(1.+beta**2.)**0.5))
    return F
