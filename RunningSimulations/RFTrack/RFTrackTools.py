import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import OctavePythonInterface as opi
import BeamDynamics as bd
try:
    import RF_Track as rft
except ModuleNotFoundError:
    pass


DEFAULT_COLOR_CYCLE = plt.rcParams["axes.prop_cycle"].by_key()['color']
DEFAULT_QUANTITIES_TO_PLOT = [
    'Bz', 'By', 'BeamPosition', 'Ez', 'mean_E', 'CaptureEfficiency',
    'Emittances', 'Sigmas', 'TwissBetas'
]


def rf_struct_from_single_period(
        fieldmapFilePath, fieldmapDim, totPeriods, powerScalingFactor, t0, phase,
        aperture=None, additionalHomogBz=None):
    rfPeriod = rf_from_field_map(
        fieldmapFilePath, fieldmapDim, powerScalingFactor, t0, phase,
        aperture=aperture, additionalHomogBz=additionalHomogBz)
    structLat = rft.Lattice()
    for rfPeriodInd in np.arange(totPeriods):
        structLat.append(rfPeriod)
    return structLat


def rf_from_field_map(
        fieldmapOrFilePath, fieldmapDim, powerScalingFactor, t0, phase, smooth=0,
        aperture=None, additionalHomogBz=None):
    try:
        fieldmapOrFilePath.keys()
        rfField = fieldmapOrFilePath
    except AttributeError:
        rfField = opi.load_octave_matrices(fieldmapOrFilePath)
    # TODO: Check generality of following line
    structL = np.diff(rfField['Z'].flatten()[[0, -1]])[0]  # [m]
    if fieldmapDim == '1D':
        dz = rfField['Z'][1] - rfField['Z'][0]  # [m]
        rf = rft.RF_FieldMap_1d(
            rfField['Ez'], dz, structL, rfField['frequency'], rfField['wave_direction'])  # [V/m]
    if fieldmapDim == '1D_3rdOrderExpansion':
        dz = rfField['Z'][1] - rfField['Z'][0]  # [m]
        rf = rft.RF_FieldMap_1d_CINT(
            rfField['Ez'], dz, structL, rfField['frequency'], rfField['wave_direction'])  # [V/m]
        rf.set_smooth(smooth)
    elif fieldmapDim == '2D':
        dr = np.diff(rfField['R'][0, [0, 1]])[0]  # [m]
        dz = np.diff(rfField['Z'][[0, 1], 0])[0]  # [m]
        rf = rft.RF_FieldMap_2d(
            rfField['Er'], rfField['Ez'], rfField['Btheta'], rfField['Bz'], dr, dz, structL,
            rfField['frequency'], rfField['wave_direction'])
    elif fieldmapDim == '3D_CylindricalSym':
        dr = rfField['R'][1] - rfField['R'][0]  # [m]
        dtheta = rfField['THETA'][1] - rfField['THETA'][0]  # [rad]
        dz = rfField['Z'][1] - rfField['Z'][0]  # [m]
        rf = rft.RF_FieldMap(
            rfField['Ex'], rfField['Ey'], rfField['Ez'],
            rfField['Bx'], rfField['By'], rfField['Bz'],
            rfField['R'][0], rfField['THETA'][0], dr, dtheta, dz, structL,
            rfField['frequency'], rfField['wave_direction']
        )
        rf.set_cylindrical(True)
    else:
        raise ValueError('Invalid fieldmapDim={:s}.'.format(fieldmapDim))
    rf.set_P_map(1.)
    if powerScalingFactor is not None:
        rf.set_P_actual(powerScalingFactor)
    if t0 is not None:
        rf.set_t0(t0)
    if phase is not None:
        rf.set_phid(phase)
    if aperture is not None:
        rf.set_aperture(aperture, aperture, 'circular')
    if additionalHomogBz is not None:
        rf.set_static_Bfield(0, 0, additionalHomogBz)
    return rf


def solenoid_from_analytical_formula(L, R_IN_COIL, R_OUT_COIL, J, extensionFactor=5., dz=1e-3):
    z = np.arange(-L*extensionFactor, L*extensionFactor, dz)
    BzOnAxis = bd.generate_solenoid_fieldmap_wilson(z, 0., R_IN_COIL, R_OUT_COIL, L/2., J)
    dz = z[1] - z[0]
    solenoid = rft.Static_Magnetic_FieldMap_1d(BzOnAxis, dz)
    return solenoid


def solenoid_from_fieldmap(fieldmapFilePath, fieldmapCurrent, setCurrent):
    solField = opi.load_octave_matrices(fieldmapFilePath)
    dz = solField['Z'][1] - solField['Z'][0]  # [m]
    BzOnAxis = solField['Bz'] / fieldmapCurrent * setCurrent  # [T]
    solenoid = rft.Static_Magnetic_FieldMap_1d(BzOnAxis, dz)
    return solenoid


def save_em_fields(
        volOrLat, xMesh, yMesh, zMesh,
        outRelPath=None, outSuffix=None, returnMultidimNpArray=False):
    emFieldsNp = np.zeros([len(zMesh)*len(yMesh)*len(xMesh), 9])
    meshInd = 0
    for x in xMesh:
        for y in yMesh:
            for z in zMesh:
                E, B = volOrLat.get_field(x*1e3, y*1e3, z*1e3, 0)
                #   x [mm], y [mm], z [mm], t [mm/c]
                emFieldsNp[meshInd, 0] = x
                emFieldsNp[meshInd, 1] = y
                emFieldsNp[meshInd, 2] = z
                emFieldsNp[meshInd, 3] = E[0]
                emFieldsNp[meshInd, 4] = E[1]
                emFieldsNp[meshInd, 5] = E[2]
                emFieldsNp[meshInd, 6] = B[0]
                emFieldsNp[meshInd, 7] = B[1]
                emFieldsNp[meshInd, 8] = B[2]
                meshInd += 1
    if returnMultidimNpArray:
        return emFieldsNp.reshape([len(xMesh), len(yMesh), len(zMesh), 9])
    emFields = pd.DataFrame(emFieldsNp, columns=['x', 'y', 'z', 'Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz'])
    emFields.to_csv(os.path.join(outRelPath, 'EMFields{:s}.dat'.format(outSuffix)), index=None)
    return emFields


def save_plot_transport(ax, volOrLat, beam0, beam1, outRelPath, outSuffix=''):
    try:
        # Volume
        # TODO: Check indexing with Andrea
        zMeshStart = volOrLat.get_s0()[0, 2] / 1e3  # [m]
        zMeshEnd = volOrLat.get_s1()[0, 2] / 1e3  # [m]
        getTransportTableStr = '%mean_S %mean_X %mean_Y %sigma_X %sigma_Y'
        sInd = 4
    except AttributeError as err:
        if err.args[0][1:8] == 'Lattice':
            zMeshStart = 0.
            zMeshEnd = volOrLat.get_length()  # [m]
            getTransportTableStr = '%S %mean_x %mean_y %sigma_x %sigma_y'
            sInd = 6
    zMesh = np.arange(zMeshStart, zMeshEnd, 1e-3)  # [m]
    # zMesh = np.linspace(zMeshStart, zMeshEnd, 1000)  # [m]
    getTransportTableStr += ' %mean_E %emitt_x %emitt_y %emitt_4d %alpha_x %beta_x %alpha_y %beta_y'
    if outSuffix != '':
        outSuffix = '_' + outSuffix
    emFields = save_em_fields(volOrLat, [0], [0], zMesh, outRelPath, outSuffix)
    # Get transport table
    TT = volOrLat.get_transport_table(getTransportTableStr)
    transportTable = pd.DataFrame(TT, columns=getTransportTableStr.replace('%', '').split())
    transportTable.to_csv(
        os.path.join(outRelPath, 'TransportTable{:s}.dat'.format(outSuffix)), index=None
    )
    # Compute capture efficiency
    M0 = beam0.get_phase_space()
    Mlost = beam1.get_lost_particles()
    # Columns of Mlost like columns 1-6 of M0, in addition:
    # t [mm/c] (Volume) or z [mm] (Lattice) at which particle was lost, m [kg],
    # Q [?] of particle type, Q of macro-particle [?]
    try:
        Mlost = Mlost[Mlost[:, sInd].argsort()]
        sCapture = Mlost[:, sInd]
        captureEff = (M0.shape[0] - np.arange(1, Mlost.shape[0]+1, 1)) / M0.shape[0]
    except IndexError:
        sCapture = np.array(TT[[0, -1], 0])
        captureEff = np.ones(sCapture.shape)
    captureEfficiency = pd.DataFrame(
        np.row_stack([sCapture, captureEff]).T, columns=['s', 'CaptureEfficiency'])
    captureEfficiency.to_csv(
        os.path.join(outRelPath, 'CaptureEfficiency{:s}.dat'.format(outSuffix)), index=None)
    plot_transport(ax, emFields, transportTable, captureEfficiency)


def load_plot_transport(
        ax, simRelPath, fileSuffix='', sShiftEMFields=0, sShiftGlobal=0, normFactorCaptureEff=1.,
        quantitiesToPlot=DEFAULT_QUANTITIES_TO_PLOT,
        lineStyles=['-']*len(DEFAULT_QUANTITIES_TO_PLOT), colors=None, cutAtInversionOfS=False):
    # TODO: Global variable, use also in save_plot_transport()
    transportFileBases = ['EMFields', 'TransportTable', 'CaptureEfficiency']
    if fileSuffix != '':
        transportFileNames = [fBase + '_' + fileSuffix + '.dat' for fBase in transportFileBases]
    else:
        transportFileNames = [fBase + '.dat' for fBase in transportFileBases]
    transportDfs = []
    for fileName in transportFileNames:
        transportDfs.append(
            pd.read_csv(os.path.join(simRelPath, fileName))
        )
    plot_transport(
        ax, *transportDfs, sShiftEMFields=sShiftEMFields, sShiftGlobal=sShiftGlobal,
        normFactorCaptureEff=normFactorCaptureEff, quantitiesToPlot=quantitiesToPlot,
        lineStyles=lineStyles, colors=colors, cutAtInversionOfS=cutAtInversionOfS)


# TODO: Move following function to module BeamDynamics?
def plot_transport(
        ax, emFields, transpTab, captureEff, sShiftEMFields=0, sShiftGlobal=0,
        normFactorCaptureEff=1., quantitiesToPlot=DEFAULT_QUANTITIES_TO_PLOT,
        lineStyles=['-']*len(DEFAULT_QUANTITIES_TO_PLOT), colors=None, cutAtInversionOfS=False):
    try:
        s = transpTab['mean_S'].to_numpy() / 1e3  # [m]
        meanXName = 'mean_X'
        meanYName = 'mean_Y'
        sigmaXName = 'sigma_X'
        sigmaYName = 'sigma_Y'
    except KeyError:
        # TODO: Verify mean_t --> S is OK for existing analyses and scripts
        s = transpTab['S'].to_numpy()  # [m]
        meanXName = 'mean_x'
        meanYName = 'mean_y'
        sigmaXName = 'sigma_x'
        sigmaYName = 'sigma_y'
    if cutAtInversionOfS:
        try:
            inversionInd = np.where(s[1:] - s[:-1] < 0)[0][0]
            s = s[:inversionInd]
            transpTab = transpTab[:inversionInd]
        except IndexError:
            pass
    # TODO: Is sShiftEMFields still necessary?
    s += sShiftGlobal
    sEmFields = emFields['z'] + sShiftEMFields + sShiftGlobal
    sCaptureEff = captureEff['s'] / 1e3 + sShiftEMFields + sShiftGlobal
    # TODO: Recognize when axis is empty (xlim = [0, 1])
    sLims = np.array([
        np.min([ax[0].get_xlim()[0], s.min(), sEmFields.min(), sCaptureEff.min()]),
        np.max([ax[0].get_xlim()[1], s.max(), sEmFields.max(), sCaptureEff.max()])])
    if lineStyles is None or isinstance(lineStyles, str):
        lineStyles = [lineStyles] * len(quantitiesToPlot)
    if colors is None or isinstance(colors, str):
        colors = [colors] * len(quantitiesToPlot)
    for axInd, (quantity, lineStyle, color) in enumerate(zip(quantitiesToPlot, lineStyles, colors)):
        if quantity == 'By':
            ax[axInd].plot(sEmFields, emFields['By'], lineStyle, color=color)
            # , color=DEFAULT_COLOR_CYCLE[0]
            ax[axInd].set_xlim(sLims)
            ByLims = np.array([
                np.min([ax[axInd].get_ylim()[0], emFields['By'].min()]),
                np.max([ax[axInd].get_ylim()[1], emFields['By'].max()])])
            ax[axInd].set_ylim(ByLims)
            ax[axInd].set_ylabel('By [T]')  # , color=DEFAULT_COLOR_CYCLE[0]
            ax[axInd].grid(True)
        elif quantity == 'Bz':
            ax[axInd].plot(sEmFields, emFields['Bz'], lineStyle, color=color)
            # , color=DEFAULT_COLOR_CYCLE[0]
            ax[axInd].set_xlim(sLims)
            BzLims = np.array([0, np.max([ax[axInd].get_ylim()[1], emFields['Bz'].max()])])
            ax[axInd].set_ylim(BzLims)
            ax[axInd].set_ylabel('Bz [T]')  # , color=DEFAULT_COLOR_CYCLE[0]
            ax[axInd].grid(True)
        elif quantity == 'Ez':
            ax[axInd].plot(sEmFields, emFields['Ez']/1e6, lineStyle, color=color)
            # , color=DEFAULT_COLOR_CYCLE[0]
            ax[axInd].set_xlim(sLims)
            EzLims = np.array([
                np.min([ax[axInd].get_ylim()[0], emFields['Ez'].min()]),
                np.max([ax[axInd].get_ylim()[1], emFields['Ez'].max()])])
            ax[axInd].set_ylim(EzLims/1e6)
            ax[axInd].set_ylabel('Ez [MV/m]')  # , color=DEFAULT_COLOR_CYCLE[0]
            ax[axInd].grid(True)
        elif quantity == 'BeamPosition':
            noMarkers = 20
            markEvery = int(len(s) / noMarkers)
            markersPerMeter = 1.  # [m]
            markEvery = int(len(s) / (s[-1] - s[0]) / markersPerMeter)
            if markEvery < 1:
                markEvery = 1
            p = ax[axInd].plot(
                s, transpTab[meanXName], lineStyle+'v', markevery=markEvery, color=color)
            ax[axInd].plot(
                s, transpTab[meanYName],
                lineStyle+'^', markevery=markEvery, color=p[0].get_color())
            ax[axInd].set_xlim(sLims)
            beamPositionLims = np.array([
                np.min([ax[axInd].get_ylim()[0], transpTab[[meanXName, meanYName]].stack().min()]),
                np.max([ax[axInd].get_ylim()[1], transpTab[[meanXName, meanYName]].stack().max()])])
            ax[axInd].set_ylim(beamPositionLims)
            ax[axInd].set_ylabel('Beam pos. [mm]')
            ax[axInd].legend(['x', 'y'])
            ax[axInd].grid(True)
        elif quantity == 'mean_E':
            ax[axInd].plot(s, transpTab['mean_E'], lineStyle, color=color)
            ax[axInd].set_xlim(sLims)
            Elims = np.array([
                np.min([ax[axInd].get_ylim()[0], transpTab['mean_E'].min()]),
                np.max([ax[axInd].get_ylim()[1], transpTab['mean_E'].max()])])
            ax[axInd].set_ylim(Elims)
            ax[axInd].set_ylabel('E [MeV]')
            ax[axInd].grid(True)
        elif quantity == 'CaptureEfficiency':
            ax[axInd].plot(
                sCaptureEff, captureEff['CaptureEfficiency']*normFactorCaptureEff,
                lineStyle, color=color)
            ax[axInd].set_xlim(sLims)
            ax[axInd].set_ylim([0, 1])
            ax[axInd].set_ylabel('Capture eff.')
            ax[axInd].grid(True)
        elif quantity == 'Emittances':
            noMarkers = 20
            markEvery = int(len(s) / noMarkers)
            markersPerMeter = 1.  # [m]
            markEvery = int(len(s) / (s[-1] - s[0]) / markersPerMeter)
            if markEvery < 1:
                markEvery = 1
            p = ax[axInd].plot(
                s, transpTab['emitt_x']/1e3, lineStyle+'v', markevery=markEvery, color=color)
            ax[axInd].plot(
                s, transpTab['emitt_y']/1e3,
                lineStyle+'^', markevery=markEvery, color=p[0].get_color())
            ax[axInd].plot(
                s, transpTab['emitt_4d']/1e3,
                lineStyle, markevery=markEvery, color=p[0].get_color())
            ax[axInd].set_xlim(sLims)
            emitLims = np.array([
                np.min([
                    ax[axInd].get_ylim()[0],
                    transpTab[['emitt_x', 'emitt_y', 'emitt_4d']].stack().min()]),
                np.max([
                    ax[axInd].get_ylim()[1],
                    transpTab[['emitt_x', 'emitt_y', 'emitt_4d']].stack().max()])])
            ax[axInd].set_ylim(emitLims/1e3)
            ax[axInd].set_ylabel('Emitt. [pimmrad]')
            ax[axInd].legend(['x (2d)', 'y (2d)', 'Trans. 4d'])
            ax[axInd].grid(True)
        elif quantity == 'Sigmas':
            noMarkers = 20
            markEvery = int(len(s) / noMarkers)
            markersPerMeter = 1.  # [m]
            markEvery = int(len(s) / (s[-1] - s[0]) / markersPerMeter)
            if markEvery < 1:
                markEvery = 1
            p = ax[axInd].plot(
                s, transpTab[sigmaXName], lineStyle+'v', markevery=markEvery, color=color)
            ax[axInd].plot(
                s, transpTab[sigmaYName],
                lineStyle+'^', markevery=markEvery, color=p[0].get_color())
            ax[axInd].set_xlim(sLims)
            sigmaLims = np.array([
                np.min([
                    ax[axInd].get_ylim()[0], transpTab[[sigmaXName, sigmaYName]].stack().min()]),
                np.max([
                    ax[axInd].get_ylim()[1], transpTab[[sigmaXName, sigmaYName]].stack().max()])])
            ax[axInd].set_ylim(sigmaLims)
            ax[axInd].set_ylabel('Sigma [mm]')
            ax[axInd].legend(['x (2d)', 'y (2d)'])
            ax[axInd].grid(True)
        elif quantity == 'TwissBetas':
            try:
                betaLims = np.array([
                    np.min([ax[axInd].get_ylim()[0], transpTab['beta_x'].min()]),
                    np.max([ax[axInd].get_ylim()[1], transpTab['beta_y'].max()])
                ])
                p = ax[axInd].plot(
                    s, transpTab['beta_x'], lineStyle+'v', markevery=markEvery, color=color)
                ax[axInd].plot(
                    s, transpTab['beta_y'],
                    lineStyle+'^', markevery=markEvery, color=p[0].get_color())
                ax[axInd].set_xlim(sLims)
                ax[axInd].set_ylim(betaLims)
                ax[axInd].set_ylabel('Betatron function [m]')
                ax[axInd].grid(True)
            except KeyError:
                pass
        elif quantity is None:
            ax[axInd].set_visible(False)
        ax[-1].set_xlabel('s [m]')
