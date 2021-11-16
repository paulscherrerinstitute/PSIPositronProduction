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
import matplotlib.markers as pltMarkers



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
UNITS_STANDARD_DF = {
    'x': 'mm', 'px': 'MeV/c',
    'y': 'mm', 'py': 'MeV/c',
    'z': 'mm', 'pz': 'MeV/c',
    't': 'ns',
    'pdgId': '', 'Q': 'C'
}
PRECISION_STANDARD_DF = 9

SCI_NOTATION_FORMATTER = ('{:' + str(PRECISION_STANDARD_DF+9) + '.' + str(PRECISION_STANDARD_DF) + 'e}').format
INTEGER_FORMATTER = ('{:'+str(PRECISION_STANDARD_DF+9)+'d}').format

FILE_TYPES_SPECS = {
    'standardDf': {
        'ext': '.sdf_txt',
        'columnOrder': COLUMN_ORDER_STANDARD_DF,
        'header': True,
        'formatters': [SCI_NOTATION_FORMATTER]*len(COLUMN_ORDER_STANDARD_DF)
            # formatterStr = '{:'+ str(PRECISION_STANDARD_DF+9) + '.' + str(PRECISION_STANDARD_DF) + 'e}'
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
        'columnOrder': ('x', 'xp', 'y', 'yp', 't', 'p', 'trackingId'),
        'header': None
    },
}

DATA_BASE_PATH = '/afs/psi.ch/project/Pcubed'


def build_data_path(relPath):
    if relPath[0] in ['/', '\\']:
        raise ValueError('relPath cannot start with a path separator ({:s} found).'.format(relPath[0]))
    absPath = os.path.join(DATA_BASE_PATH, relPath)
    return absPath


def pdgId_to_particle_const(pdgId, constName):
    flagSingleValue = False
    if not isinstance(pdgId, (pd.Series,np.ndarray)):
        flagSingleValue = True
        pdgId = [pdgId]
    Erest = pd.Series([
        PART_CONSTS[constName][pid] if pid in PART_CONSTS[constName].keys() else np.nan for pid in pdgId
    ])
    if flagSingleValue:
        Erest = Erest[0]
    return Erest


def p_to_beta(p, pdgId):
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    beta = np.sqrt(p**2. / (p**2. + Erest**2.))
    return beta


def pVect_to_p(px, py, pz):
    p = np.sqrt(px**2. + py**2. + pz**2.)
    return p


def pTransv_to_slope(pTransv, pLong):
    slope = np.arctan(pTransv / pLong) * 1e3                                    # [mrad]
    return slope


def p_to_Ekin(p, pdgId):
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    E = np.sqrt(p**2. + Erest**2.)
    Ekin = E - Erest
    return Ekin


def Ekin_to_p(Ekin, pdgId):
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    p = np.sqrt((Ekin+Erest)**2. - Erest**2.)
    return p


def Ekin_to_gamma(Ekin, pdgId):
    Erest = pdgId_to_particle_const(pdgId, 'Erest')
    gamma = Ekin / Erest + 1
    return gamma


def gamma_to_beta(gamma):
    beta = np.sqrt(1. - (1./gamma**2.))
    return beta


def Ekin_to_beta(Ekin, pdgId):
    gamma = Ekin_to_gamma(Ekin, pdgId)
    beta = gamma_to_beta(gamma)
    return beta


def extend_standard_df(standardDf):
    p = pVect_to_p(standardDf['px'], standardDf['py'], standardDf['pz'])
    dfExtension = pd.DataFrame()
    dfExtension['Ekin'] = p_to_Ekin(p, standardDf['pdgId'])
    dfExtension['gammaRel'] = Ekin_to_gamma(dfExtension['Ekin'], standardDf['pdgId'])
    dfExtension['betaRel'] = p_to_beta(p, standardDf['pdgId'])
    dfExtension['xp'] = pTransv_to_slope(standardDf['px'], standardDf['pz'])
    dfExtension['yp'] = pTransv_to_slope(standardDf['py'], standardDf['pz'])
    for colName in dfExtension.columns:
        if colName in standardDf.columns:
            warnings.warn('Overwriting available column {:s}.'.format(colName))
            relDiffs = (standardDf[colName] - dfExtension[colName]) / dfExtension[colName]
            relDiffThreshold = 1e-6
            isDiff = np.abs(relDiffs) > relDiffThreshold
            if np.any(isDiff):
                warnings.warn(
                    '{:d} discrepancies larger than {:e} detected between available and recomputed values.'.format(
                        np.sum(isDiff), relDiffThreshold
                ))
    standardDf = dfExtension.combine_first(standardDf)
    return standardDf


def compute_emittance(standardDf, planeName, correctOffsets=False):
    x = standardDf[planeName]
    px = standardDf['p'+planeName]
    thresholdFactor = 0.001
    xAvg = x.mean()
    if xAvg > x.std()*thresholdFactor:
        warnings.warn(
            'Average position {:s}Avg = {:.3f} mm.'.format(planeName, xAvg)
        )
    pxAvg = px.mean()
    if pxAvg > px.std()*thresholdFactor:
        warnings.warn(
            'Average divergence p{:s}Avg = {:.3f} MeV/c.'.format(
                planeName, pxAvg
        ))
    if correctOffsets:
        x -= xAvg
        px -= pxAvg
        print(
            'Correcting offsets {0:s}Avg = {1:.3f} mm and p{0:s}Avg = {2:.3f} MeV/c.'.format(
                planeName, xAvg, pxAvg
        ))
    # TODO: Check that all pdgIds are identical
    Erest = pdgId_to_particle_const(standardDf['pdgId'][0], 'Erest')
    emitNorm = np.sqrt(
        (x**2.).sum() * (px**2.).sum() - (x*px).sum()**2.
    ) / x.shape[0] / Erest * 1e3
    return emitNorm                                                             # [mm mrad]


def export_sim_db_to_xlsx(jsonFilePath):
    dbFile = open(jsonFilePath)
    dbJson = json.load(dbFile)
    dbDf = pd.json_normalize(dbJson)
    xlsxFilePath = os.path.splitext(jsonFilePath)[0] + '.xlsx'
    dbDf.to_excel(xlsxFilePath)


# TODO: Changed from sourceFilePath to outFilePath=None
    # outFilePath = sourceFilePath + additionalLabel + fileTypeSpecs['ext']
def generate_fwf(df, formatType='standardDf', outFilePath=None):
    fileTypeSpecs = FILE_TYPES_SPECS[formatType]
    dfStr = df[fileTypeSpecs['columnOrder']].to_string(
        formatters=fileTypeSpecs['formatters'], index=False,
        header=fileTypeSpecs['header']
    )
    # headerList = standardDf.columns + ' [' + UNITS_STANDARD_DF + ']'
    if outFilePath is not None:
        outFilePath += FILE_TYPES_SPECS[formatType]['ext']
        with open(outFilePath, 'w') as outFile:
            outFile.write(dfStr)
    return dfStr


def load_standard_fwf(sourceFilePath):
    standardDf = pd.read_fwf(sourceFilePath)
    standardDf = extend_standard_df(standardDf)
    return standardDf


def convert_irina_distr_to_standard_df(sourceFilePath, saveStandardFwf=False):
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
    if saveStandardFwf:
        generate_fwf(standardDf, outFilePath=sourceFilePath)
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
        standardDf['t'] = standardDf['t'] * 1e-3                                # [ns]
        standardDf['pdgId'] = standardDf['pdgId'].astype(int)
        standardDf['Q'] = pdgId_to_particle_const(standardDf['pdgId'], 'Q')
        standardDf = extend_standard_df(standardDf)
        # TODO: Check if 'e' in FCC-ee Target Tracking is Ekin or Ekin + Erest (simple curiosity)
        # One could have also used: standardDf.rename(columns={'e': 'Ekin'}, inplace=True)
        fileSuffix = '_' + distrName
        if distrName == 'amor_leave' and pdgId:
            standardDf = standardDf[standardDf['pdgId'].isin(pdgId)]
            fileSuffix += '_pdgId'
            for id in pdgId:
                fileSuffix += '_' + str(id)
        if saveStandardFwf:
            filePath, fileExt = os.path.splitext(sourceFilePath)
            generate_fwf(standardDf, outFilePath=filePath+fileSuffix+fileExt)
        dfDict[distrName] = standardDf
    return dfDict


def convert_astra_to_standard_df(sourceFilePath, zProjection=None, zCut=None, saveStandardFwf=False, verbose=False):
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
    standardDf.loc[1:,'z'] = standardDf.loc[1:,'z'] + standardDf.loc[0,'z']
    standardDf.loc[1:,'pz'] = standardDf.loc[1:,'pz'] + standardDf.loc[0,'pz']
    standardDf.loc[1:,'clock'] = \
        standardDf.loc[1:,'clock'] + standardDf.loc[0,'clock']                  # [ns]
    standardDf.rename(columns={'clock': 't'}, inplace=True)
    standardDf['particleIndex'].replace(to_replace={1:11, 2:-11}, inplace=True)
    standardDf.rename(columns={'particleIndex': 'pdgId'}, inplace=True)
    standardDf['macroCharge'] = standardDf['macroCharge'] * 1.e-9               # [C]
    standardDf.rename(columns={'macroCharge': 'Q'}, inplace=True)
    standardDf.drop('statusFlag', axis=1, inplace=True)
    standardDf = extend_standard_df(standardDf)
    labelStr = ''
    # If requested, cut particles that did not reach zCut
    if zCut is not None:
        cutInds = standardDf['z'] < zCut
        rejectedParticles = standardDf[cutInds]
        standardDf = standardDf[~cutInds]
        warnings.warn(
            '{:d} particles with z < {:.3f} mm have been discarded.'.format(
                cutInds.sum(), zCut
        ))
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
    if saveStandardFwf:
        generate_fwf(standardDf, outFilePath=sourceFilePath+labelStr)
    return standardDf


def convert_standard_df_to_astra(standardDf=None, sourceFilePath=None, refParticleId=0, saveAstraDistr=False, outFilePath='AstraDistribution'):
    # TODO: saveAstraDistr and outFilePath redundant, check also similar functions
    standardDf, outFilePath = convert_from_standard_df_input_check(
        standardDf, sourceFilePath, outFilePath
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
    astraDf['pdgId'].replace(to_replace={11:1, -11:2}, inplace=True)
    # Small check for the reference particle
    if not np.isclose(astraDf['x'][refParticleId], 0) or not np.isclose(astraDf['y'][refParticleId], 0):
        warnings.warn('Reference particle is out of axis with (x,y) = ({:.6f}, {:.6f}) m'.format(
            astraDf['x'][refParticleId], astraDf['y'][refParticleId]
    ))
    # Longitudinal variables with respect to the reference particle
    astraDf.loc[1:,'z'] = astraDf.loc[1:,'z'] - astraDf.loc[refParticleId,'z']
    astraDf.loc[1:,'pz'] = astraDf.loc[1:,'pz'] - astraDf.loc[refParticleId,'pz']
    astraDf.loc[1:,'t'] = astraDf.loc[1:,'t'] - astraDf.loc[refParticleId,'t']
    # Add status flag
    statusFlag = np.full(astraDf.shape[0], 5)
    astraDf['statusFlag'] = statusFlag
    astraDf['pdgId'] = astraDf['pdgId'].astype(int)
    astraDf['statusFlag'] = astraDf['statusFlag'].astype(int)
    # TODO: Move reference particle to first position in the data frame
    astraDf.columns = FILE_TYPES_SPECS['astra']['columnOrder']
    if saveAstraDistr:
        generate_fwf(astraDf, formatType='astra', outFilePath=outFilePath)
    return astraDf


def convert_from_standard_df_input_check(standardDf, sourceFilePath, outFilePath):
    if sourceFilePath is not None:
        if standardDf is None:
            standardDf = load_standard_fwf(sourceFilePath)
            outFilePath = sourceFilePath
        else:
            raise ValueError(
                'Only one between standardDf and sourceFilePath can be set different than None.'
            )
    return standardDf, outFilePath


def convert_sdds_to_standard_df(sourceFilePath, z0=0, pdgId=-11, Qbunch=np.nan, saveStandardFwf=False):
    os.system('sddsconvert -ascii ' + sourceFilePath)
    standardDf = pd.read_csv(
        sourceFilePath, skiprows=24, delim_whitespace=True,
        names=FILE_TYPES_SPECS['sdds']['columnOrder']
    )
    standardDf.drop('trackingId', axis=1, inplace=True)
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
    standardDf['pdgId'] = pdgId
    Erest = pdgId_to_particle_const(standardDf['pdgId'], 'Erest')
    p = standardDf['p'] * Erest                                                 # [MeV/c]
    standardDf.drop('p', axis=1, inplace=True)
    standardDf['pz'] = p / np.sqrt(1 + np.tan(standardDf['xp'])**2. + np.tan(standardDf['yp'])**2.) # [MeV/c]
    standardDf['x'] = standardDf['x'] * 1.e3                                    # [mm]
    standardDf['px'] = np.tan(standardDf['xp']) * standardDf['pz']              # [MeV/c]
    standardDf['y'] = standardDf['y'] * 1.e3                                    # [mm]
    standardDf['py'] = np.tan(standardDf['yp']) * standardDf['pz']              # [MeV/c]
    standardDf['z'] = z0                                                        # [mm]
    standardDf['t'] = standardDf['t'] * 1.e9                                    # [ns]
    standardDf['Q'] = Qbunch / NparticlesPerBunch                               # [C]
    # Extended variables already available
    standardDf['xp'] = standardDf['xp'] * 1.e3                                  # [mrad]
    standardDf['yp'] = standardDf['yp'] * 1.e3                                  # [mrad]
    standardDf = extend_standard_df(standardDf)
    # sElegant = standardDf['betaRel'] * C * standardDf['t'] * 1e-6             # [mm]
    if saveStandardFwf:
        generate_fwf(standardDf, outFilePath=sourceFilePath)
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


def convert_standard_df_to_sdds(standardDf=None, sourceFilePath=None, outFilePath='SddsDistribution'):
    # standardDf, outFilePath = convert_from_standard_df_input_check(
    #     standardDf, sourceFilePath
    # )
    astraDf = convert_standard_df_to_astra(standardDf, sourceFilePath)
    astraDfStr = generate_fwf(astraDf, formatType='astra', outFilePath='TmpAstra')
    outFilePath += FILE_TYPES_SPECS['sdds']['ext']
    os.system('astra2elegant TmpAstra.001 ' + outFilePath)
    os.system('rm TmpAstra.001')
    # os.system('astra2elegant -pipe=in ' + astraDfStr + ' ' + outFilePath)


def plot_hist(ax, distr, binWidth=None, binLims=None, legendLabel='', orientation='vertical', parsInLabel=True, opacityHist=1.):
    defaultBinNum = 100
    if binLims is None:
        if binWidth is None:
            binWidth = (np.max(distr) - np.min(distr)) / (defaultBinNum-1)
            if binWidth == 0:
                binWidth = 1.
                defaultBinNum = 3
        binLims = [np.min(distr)-binWidth*1.5, np.max(distr)+binWidth*1.5]
    if binWidth is None:
        binEdges = np.linspace(binLims[0], binLims[1], num=defaultBinNum)
    else:
        binEdges = np.arange(binLims[0], binLims[1], binWidth)
    counts, _, histObj = ax.hist(
        distr, bins=binEdges, orientation=orientation, alpha=opacityHist
    )
    avg = np.mean(distr)
    std = np.std(distr)
    if parsInLabel:
        if legendLabel != '':
            legendLabel += ', '
        legendLabel += 'avg = {:.3f}, std = {:.3f}'.format(avg, std)
    histObj.set_label(legendLabel)
    ax.legend()
    return avg, std


def set_lims(ax, direction, var, lims):
    if lims is None:
        lims = (np.min(var), np.max(var))
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
        markerObj = pltMarkers.MarkerStyle(m)
        markerPath = markerObj.get_path().transformed(markerObj.get_transform())
        markerPaths.append(markerPath)
    pathCollectionObj.set_paths(markerPaths)


def plot_phase_space_2d(
        ax, distr, varName1=None, varName2=None, title=None, legendLabel='',
        binWidth1=None, binWidth2=None, lims1=None, lims2=None, pzCutoff=None,
        markerStyle=None, markerSize=15, color='b', opacityHist=1.
):
    if varName1 is None or varName2 is None:
        raise ValueError('varName1 and varName2 are mandatory arguments.')
    opacityScatter = 1.
    pathColl = ax[0,0].scatter(
        distr[varName1], distr[varName2],
        s=markerSize, edgecolors=color,
        facecolors='none', linewidths=1., alpha=opacityScatter
    )
    scatter_individual_marker_style(pathColl, markerStyle)
    plot_hist(
        ax[1,0], distr[varName1], binWidth=binWidth1, binLims=lims1,
        legendLabel=legendLabel, opacityHist=opacityHist
    )
    plot_hist(
        ax[0,1], distr[varName2], binWidth=binWidth2, binLims=lims2,
        legendLabel=legendLabel, orientation='horizontal',
        opacityHist=opacityHist
    )
    pzLabels = ['All, Counts = {:d}'.format(distr.shape[0]), ]
    if pzCutoff is not None:
        distrLowPz = distr[distr['pz']<pzCutoff]
        ax[0,0].scatter(
            distrLowPz[varName1], distrLowPz[varName2],
            s=markerSize, marker=markerStyle, facecolors='none', linewidths=1., edgecolors=color,
            alpha=opacityScatter
        )
        plot_hist(
            ax[1,0], distrLowPz[varName1],
            binWidth=binWidth1, binLims=lims1, opacityHist=opacityHist
        )
        plot_hist(
            ax[0,1], distrLowPz[varName2], orientation='horizontal',
            binWidth=binWidth2, binLims=lims2, opacityHist=opacityHist
        )
        pzLabels += ['pz <= {:.1f} {:s}, Counts = {:d}'.format(
            pzCutoff, UNITS_STANDARD_DF['pz'], distrLowPz.shape[0]
        )]
        ax[0,0].legend(
            pzLabels, markerscale=5.,
            loc='upper left', bbox_to_anchor=(1.2, -0.2)
        )
    set_lims(ax[0,0], 'x', distr[varName1], lims1)
    set_lims(ax[0,0], 'y', distr[varName2], lims2)
    ax[1,0].set_xlim(ax[0,0].get_xlim())
    ax[1,0].invert_yaxis()
    ax[0,1].set_ylim(ax[0,0].get_ylim())
    ax[0,0].grid(True)
    ax[0,1].grid(True)
    ax[1,0].grid(True)
    ax[0,0].set_xlabel(varName1+' ['+UNITS_STANDARD_DF[varName1]+']')
    ax[0,0].set_ylabel(varName2+' ['+UNITS_STANDARD_DF[varName2]+']')
    ax[0,0].yaxis.set_label_position('right')
    ax[0,0].yaxis.tick_right()
    ax[1,0].set_ylabel('Bin counts')
    ax[0,1].set_xlabel('Bin counts')
    if title is not None:
        ax[0,0].get_figure().suptitle(title)
    ax[1,1].set_visible(False)


def plot_distr(
        distributions, plotDefs, markerStyle=None, markerSize=15,
        title=None, legendLabels=None, figHeight=9, figWidth=16
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
    # TODO: Colors
    defaultColors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    colors = defaultColors[:len(distributions)]
    axList = []
    for plotDef in plotDefs:
        fig, ax = plt.subplots(2, 2)
        fig.set_figheight(figHeight)
        fig.set_figwidth(figWidth)
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
        # plt.show()
    return axList


def check_marker_specs(distributions, markerSpecs, specName):
    try:
        if len(markerSpecs) != len(distributions):
            raise ValueError(
                '{:s} and distributions must have the same length.'.format(
                    specName
            ))
    except TypeError:
        markerSpecs = [markerSpecs] * len(distributions)
    return markerSpecs


def generate_fieldmap_astra_ideal_tw(fileBasePath, freq, Lstructure, zRes):
    z = np.arange(0, Lstructure+zRes, zRes)
    # k = 2*np.pi*freq / C
    # Ereal = np.cos(k*z)
    # Eimag = np.sin(k*z)
    # np.savetxt(fileBasePath+'_Real.astra', np.stack([z, Ereal], axis=1), fmt='%.9e')
    # np.savetxt(fileBasePath+'_Imag.astra', np.stack([z, Eimag], axis=1), fmt='%.9e')
    Eabs = np.array([1.0]*len(z))
    np.savetxt(fileBasePath+'.astra', np.stack([z, Eabs], axis=1), fmt='%.9e')
    
    
def generate_lattice_quad_over_rf_elegant(elementName, L, Nslices, quadGradient, rfFreq, rfPhase, rfVoltage, EkinIni, pdgId=-11, ZoverA=1., elegantInputFilePath=None, quadOrder=2):
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
        inputStr += '{0:s}.Slice{1:d}: Line=({0:s}.Rf.HalfSlice, {0:s}.Nd.HalfSlice, {0:s}.Quad.{1:d}, {0:s}.Nd.HalfSlice, {0:s}.Rf.HalfSlice)\n'.format(
            elementName, sliceInd
        )
        # inputStr += generate_quad_over_rf_slice(Lslice, sliceInd, quadStrengthSlice, rfFreq, rfPhase, rfVoltageSlice)
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


def quad_strength(quadGradient, p0, pdgId=-11, ZoverA=1.):
    beta = p_to_beta(p0, pdgId)
    # TODO: Check this formula, probably without beta?
    quadStrength = C * 1e-6 * ZoverA * quadGradient / (
        beta * p0
    )
    return quadStrength


def quad_matrix(k, l):
    sqrtK = np.sqrt(np.abs(k))
    phi = sqrtK * l
    MquadFocusing = np.array((
        (np.cos(phi), np.sin(phi)/sqrtK),
        (-np.sin(phi)*sqrtK, np.cos(phi))
    ))
    MquadDefocusing = np.array((
        (np.cosh(phi), np.sinh(phi)/sqrtK),
        (np.sinh(phi)*sqrtK, np.cosh(phi))
    ))
    return MquadFocusing, MquadDefocusing


def generate_cross_distribution(xMax, yMax, p0, pzDelta, xPoints=5, yPoints=5, pzPoints=5, saveStandardFwf=False, outFilePath='CrossDistribution'):
    if not isinstance(xPoints,int) or xPoints < 1:
        raise ValueError('xPoints must be an integer >= 1.')
    if not isinstance(yPoints,int) or yPoints < 1:
        raise ValueError('yPoints must be an integer >= 1.')
    if not isinstance(pzPoints,int) or pzPoints < 1:
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
    standardDf = extend_standard_df(standardDf)
    # TODO: Remove following lines when working (this was the first implementation)
    # p = pVect_to_p(pxArray, pyArray, pzArray)
    # EkinArray = p_to_Ekin(p, pdgId)
    # gammaRelArray = Ekin_to_gamma(EkinArray, pdgId)
    # betaRelArray = Ekin_to_beta(EkinArray, pdgId)
    # xpArray = pxArray / pzArray
    # ypArray = pyArray / pzArray
    if saveStandardFwf:
        generate_fwf(standardDf, outFilePath=outFilePath)
    return standardDf
    