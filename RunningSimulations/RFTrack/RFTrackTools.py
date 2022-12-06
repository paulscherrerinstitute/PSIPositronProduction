import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import BeamDynamics as bd
try:
    import RF_Track as rft
except ModuleNotFoundError:
    pass


DEFAULT_COLOR_CYCLE = plt.rcParams["axes.prop_cycle"].by_key()['color']
DEFAULT_QUANTITIES_TO_PLOT = [
    'Bz', 'Ez', 'mean_E', 'CaptureEfficiency', 'Emittances', 'Sigmas', 'TwissBetas'
]

# TODO: Make path relative
RF_CLIC_FIELDMAP_BASEPATH = \
    '/home/tia/Repos/GIT_PSIPositronProduction/' + \
    'RFTrack/YongkeTool_V1/field/field_map_CLIC_Lband'
RF_CLIC_FREQ = 1.9986163867e+09   # [Hz]
RF_CLIC_GRADIENT = 11.23e6   # [V/m]
RF_CLIC_CELLS_PER_PERIOD = 3.
RF_CLIC_L_CELL = bd.C / RF_CLIC_FREQ / RF_CLIC_CELLS_PER_PERIOD  # [m]


def load_fieldmap_rf_clic():
    # TODO: Integrate reading of a 3 dimensional array in load_octave_matrices()
    # and thenuse that function
    meshAxes = ['ra', 'ta', 'za']
    meshes = {}
    for meshAx in meshAxes:
        filePath = RF_CLIC_FIELDMAP_BASEPATH + '_' + meshAx + '.dat'
        meshes[meshAx] = np.loadtxt(filePath, skiprows=5)
    matDimensions = [meshes[meshAx].shape[0] for meshAx in reversed(meshAxes)]
    complexComponents = ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz']
    rfFields = {}
    for complComp in complexComponents:
        filePath = RF_CLIC_FIELDMAP_BASEPATH + '_' + complComp + '.dat'
        rfFields[complComp] = np.loadtxt(
            filePath, skiprows=5, dtype=np.complex128,
            converters={0: lambda s: np.fromstring(
                s.decode("latin1").replace('(', '').replace(')', ''), sep=','
            ).view(np.complex128)}
        )
        rfFields[complComp] = rfFields[complComp].reshape(matDimensions).transpose()
        # .astype(np.complex128)
    return rfFields, meshes


def rf_clic_single_period(rfFieldmapDim):
    # Fields of 1 RF period (3 cells)
    # TODO: Data location and path construction
    rfFields, rfMesh = load_fieldmap_rf_clic()
    if rfFieldmapDim == '1D':
        rfPeriod = rft.RF_FieldMap_1d(
            rfFields['Ez'][0, 0, :],
            rfMesh['za'][1] - rfMesh['za'][0],
            rfMesh['za'][-1],
            RF_CLIC_FREQ, +1
        )
    elif rfFieldmapDim == '3D':
        rfPeriod = rft.RF_FieldMap(
            rfFields['Ex'], rfFields['Ey'], rfFields['Ez'],
            rfFields['Bx'], rfFields['By'], rfFields['Bz'],
            rfMesh['ra'][0], rfMesh['ta'][0],
            rfMesh['ra'][1] - rfMesh['ra'][0],
            rfMesh['ta'][1] - rfMesh['ta'][0],
            rfMesh['za'][1] - rfMesh['za'][0],
            rfMesh['za'][-1],
            RF_CLIC_FREQ, +1
        )
        rfPeriod.set_cylindrical(True)
    else:
        raise ValueError('Invalid rfFieldmapDim={:s}.'.format(rfFieldmapDim))
    return rfPeriod


def rf_struct_from_single_period(
        fieldmapFilePath, fieldmapDim, totPeriods, powerScalingFactor, t0, phase,
        aperture=None, additionalHomogBz=None
):
    rfPeriod = rf_clic_single_period(fieldmapDim)
    rfPeriod.set_P_map(1.)
    rfPeriod.set_P_actual(powerScalingFactor)
    if t0 is not None:
        rfPeriod.set_t0(t0)
    rfPeriod.set_phid(phase)
    if aperture is not None:
        rfPeriod.set_aperture(aperture, aperture, 'circular')
    if additionalHomogBz is not None:
        rfPeriod.set_static_Bfield(0, 0, additionalHomogBz)
    structLat = rft.Lattice()
    for rfPeriodInd in np.arange(totPeriods):
        structLat.append(rfPeriod)
    return structLat


def rf_struct_from_full_fieldmap(
        fieldmapFilePath, fieldmapDim, powerScalingFactor, t0, phase,
        aperture=None, additionalHomogBz=None
):
    rfField = bd.load_octave_matrices(fieldmapFilePath)
    dz = (rfField['Z'][1] - rfField['Z'][0]) * 1e-3  # [m]
    structL = (rfField['Z'][-1] - rfField['Z'][0]) * 1e-3  # [m]
    rfStruct = rft.RF_FieldMap_1d(rfField['E'], dz, structL, rfField['frequency'], +1)
    rfStruct.set_P_map(1.)
    if powerScalingFactor is not None:
        rfStruct.set_P_actual(powerScalingFactor)
    if t0 is not None:
        rfStruct.set_t0(t0)
    if phase is not None:
        rfStruct.set_phid(phase)
    if aperture is not None:
        rfStruct.set_aperture(aperture, aperture, 'circular')
    if additionalHomogBz is not None:
        rfStruct.set_static_Bfield(0, 0, additionalHomogBz)
    return rfStruct


def solenoid_from_analytical_formula(L, R_IN_COIL, R_OUT_COIL, J, extensionFactor=5., dz=1e-3):
    z = np.arange(-L*extensionFactor, L*extensionFactor, dz)
    BzOnAxis = bd.generate_solenoid_fieldmap_wilson(z, 0., R_IN_COIL, R_OUT_COIL, L/2., J)
    dz = z[1] - z[0]
    solenoid = rft.Static_Magnetic_FieldMap_1d(BzOnAxis, dz)
    return solenoid


def solenoid_from_fieldmap(fieldmapFilePath, fieldmapCurrent, setCurrent):
    solField = bd.load_octave_matrices(fieldmapFilePath)
    dz = solField['Z'][1] - solField['Z'][0]  # [m]
    BzOnAxis = solField['Bz'] / fieldmapCurrent * setCurrent  # [T]
    solenoid = rft.Static_Magnetic_FieldMap_1d(BzOnAxis, dz)
    return solenoid


def save_plot_transport(ax, vol, beam0, beam1, outRelPath, outSuffix=''):
    # Prepare Ez and Bz for plotting
    zAxis = np.arange(vol.get_s0(), vol.get_s1(), 1e-3)   # [m]
    # zAxis = np.linspace(vol.get_s0(), vol.get_s1(), 1000)   # [m]
    Ez = []
    Bz = []
    for z in zAxis:
        E, B = vol.get_field(0, 0, z*1e3, 0)   # x,y,z,t (mm, mm/c)
        Ez.append(E[2])
        Bz.append(B[2])
    Ez = np.array(Ez)
    Bz = np.array(Bz)
    emFields = pd.DataFrame(np.row_stack([zAxis, Ez, Bz]).T, columns=['z', 'Ez', 'Bz'])
    if outSuffix != '':
        outSuffix = '_' + outSuffix
    emFields.to_csv(os.path.join(outRelPath, 'EMFields{:s}.dat'.format(outSuffix)), index=None)
    # Get transport table
    getTransportTableStr = '%mean_S %emitt_x %emitt_y %emitt_4d %sigma_X %sigma_Y %mean_E'
    TT = vol.get_transport_table(getTransportTableStr)
    transportTable = pd.DataFrame(TT, columns=getTransportTableStr.replace('%', '').split())
    transportTable.to_csv(
        os.path.join(outRelPath, 'TransportTable{:s}.dat'.format(outSuffix)), index=None
    )
    # Compute capture efficiency
    M0 = beam0.get_phase_space()
    Mlost = beam1.get_lost_particles()
    # Columns of Mlost like columns 1-6 of M0, in addition:
    # t [mm/c] at which particle was lost, m [kg],
    # Q [?] of particle type, Q of macro-particle [?]
    try:
        Mlost = Mlost[Mlost[:, 4].argsort()]
        sCapture = Mlost[:, 4]
        captureEff = (
            M0.shape[0] - np.arange(1, Mlost.shape[0]+1, 1)
        ) / M0.shape[0]
    except IndexError:
        sCapture = np.array(TT[[0, -1], 0])
        captureEff = np.ones(sCapture.shape)
    captureEfficiency = pd.DataFrame(
        np.row_stack([sCapture, captureEff]).T, columns=['s', 'CaptureEfficiency']
    )
    captureEfficiency.to_csv(
        os.path.join(outRelPath, 'CaptureEfficiency{:s}.dat'.format(outSuffix)), index=None
    )
    plot_transport(ax, emFields, transportTable, captureEfficiency)


def save_plot_transport_in_lattice(ax, lat, beam0, beam1, outRelPath, outSuffix=''):
    # TODO: Integrate in save_plot_transport for volume, differences are small.
    # Prepare Ez and Bz for plotting
    zAxis = np.arange(0., lat.get_length(), 1e-3)   # [m]
    Ez = []
    Bz = []
    for z in zAxis:
        E, B = lat.get_field(0, 0, z*1e3, 0)   # x,y,z,t (mm, mm/c)
        Ez.append(E[2])
        Bz.append(B[2])
    Ez = np.array(Ez)
    Bz = np.array(Bz)
    emFields = pd.DataFrame(np.row_stack([zAxis, Ez, Bz]).T, columns=['z', 'Ez', 'Bz'])
    if outSuffix != '':
        outSuffix = '_' + outSuffix
    emFields.to_csv(os.path.join(outRelPath, 'EMFields{:s}.dat'.format(outSuffix)), index=None)
    # Get transport table
    # TODO: Add Twiss parameters
    getTransportTableStr = '%S %emitt_x %emitt_y %emitt_4d %sigma_x %sigma_y %mean_E ' + \
        '%beta_x %beta_y'
    TT = lat.get_transport_table(getTransportTableStr)
    transportTable = pd.DataFrame(TT, columns=getTransportTableStr.replace('%', '').split())
    transportTable.to_csv(
        os.path.join(outRelPath, 'TransportTable{:s}.dat'.format(outSuffix)), index=None
    )
    # Compute capture efficiency
    M0 = beam0.get_phase_space()
    Mlost = beam1.get_lost_particles()
    # Columns of Mlost like columns 1-6 of M0, in addition:
    # t [mm/c] at which particle was lost, m [kg],
    # Q [?] of particle type, Q of macro-particle [?]
    try:
        # sInd = 4  # in Volume()
        sInd = 6  # in Lattice()
        Mlost = Mlost[Mlost[:, sInd].argsort()]
        sCapture = Mlost[:, sInd]
        captureEff = (
            M0.shape[0] - np.arange(1, Mlost.shape[0]+1, 1)
        ) / M0.shape[0]
    except IndexError:
        sCapture = np.array(TT[[0, -1], 0])
        captureEff = np.ones(sCapture.shape)
    captureEfficiency = pd.DataFrame(
        np.row_stack([sCapture, captureEff]).T, columns=['s', 'CaptureEfficiency']
    )
    captureEfficiency.to_csv(
        os.path.join(outRelPath, 'CaptureEfficiency{:s}.dat'.format(outSuffix)), index=None
    )
    plot_transport(ax, emFields, transportTable, captureEfficiency)


def load_plot_transport(
        ax, simRelPath, fileSuffix='', sShiftEMFields=0, sShiftGlobal=0, normFactorCaptureEff=1.,
        quantitiesToPlot=DEFAULT_QUANTITIES_TO_PLOT
):
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
        normFactorCaptureEff=normFactorCaptureEff, quantitiesToPlot=quantitiesToPlot)


# TODO: Move following function to module BeamDynamics?
def plot_transport(
        ax, emFields, transpTab, captureEff, sShiftEMFields=0, sShiftGlobal=0,
        normFactorCaptureEff=1., quantitiesToPlot=DEFAULT_QUANTITIES_TO_PLOT):
    try:
        s = transpTab['mean_S'] / 1e3  # [m]
        sigmaXName = 'sigma_X'
        sigmaYName = 'sigma_Y'
    except KeyError:
        # TODO: Verify mean_t --> S is OK for existing analyses and scripts
        s = transpTab['S']  # [m]
        sigmaXName = 'sigma_x'
        sigmaYName = 'sigma_y'
    # TODO: Is sShiftEMFields still necessary?
    s += sShiftGlobal
    # TODO: Recognize when axis is empty (xlim = [0, 1])
    sLims = np.array([
        np.min([ax[0].get_xlim()[0], s.min()]),
        np.max([ax[0].get_xlim()[1], s.max()])
    ])
    for axInd, quantity in enumerate(quantitiesToPlot):
        if quantity == 'Bz':
            ax[axInd].plot(emFields['z']+sShiftEMFields+sShiftGlobal, emFields['Bz'])
            # , color=DEFAULT_COLOR_CYCLE[0]
            ax[axInd].set_xlim(sLims)
            BzLims = np.array([0, np.max([ax[axInd].get_ylim()[1], emFields['Bz'].max()])])
            ax[axInd].set_ylim(BzLims)
            ax[axInd].set_xlabel('s [m]')
            ax[axInd].set_ylabel('Bz [T]')  # , color=DEFAULT_COLOR_CYCLE[0]
            ax[axInd].grid(True)
        elif quantity == 'Ez':
            ax[axInd].plot(emFields['z']+sShiftEMFields+sShiftGlobal, emFields['Ez']/1e6, '--')
            # , color=DEFAULT_COLOR_CYCLE[0]
            ax[axInd].set_xlim(sLims)
            EzLims = np.array([
                np.min([ax[axInd].get_ylim()[0], emFields['Ez'].min()]),
                np.max([ax[axInd].get_ylim()[1], emFields['Ez'].max()])])
            ax[axInd].set_ylim(EzLims/1e6)
            ax[axInd].set_ylabel('Ez [MV/m]')  # , color=DEFAULT_COLOR_CYCLE[0]
            ax[axInd].grid(True)
        elif quantity == 'mean_E':
            ax[axInd].plot(s, transpTab['mean_E'])
            ax[axInd].set_xlim(sLims)
            Elims = np.array([
                np.min([ax[axInd].get_ylim()[0], transpTab['mean_E'].min()]),
                np.max([ax[axInd].get_ylim()[1], transpTab['mean_E'].max()])])
            ax[axInd].set_ylim(Elims)
            ax[axInd].set_xlabel('s [m]')
            ax[axInd].set_ylabel('E [MeV]')
            ax[axInd].grid(True)
        elif quantity == 'CaptureEfficiency':
            ax[axInd].plot(
                captureEff['s']/1e3+sShiftEMFields+sShiftGlobal,
                captureEff['CaptureEfficiency']*normFactorCaptureEff, '--'
            )
            ax[axInd].set_xlim(sLims)
            ax[axInd].set_ylim([0, 1])
            ax[axInd].set_ylabel('Capture eff.')
            ax[axInd].grid(True)
        elif quantity == 'Emittances':
            noMarkers = 20
            markEvery = int(len(s) / noMarkers)
            if markEvery < 1:
                markEvery = 1
            p = ax[axInd].plot(s, transpTab['emitt_x']/1e3, '-v', markevery=markEvery)
            ax[axInd].plot(
                s, transpTab['emitt_y']/1e3, '-^', markevery=markEvery, color=p[0].get_color())
            ax[axInd].plot(
                s, transpTab['emitt_4d']/1e3, '-', markevery=markEvery, color=p[0].get_color())
            ax[axInd].set_xlim(sLims)
            emitLims = np.array([
                np.min([
                    ax[axInd].get_ylim()[0],
                    transpTab[['emitt_x', 'emitt_y', 'emitt_4d']].stack().min()]),
                np.max([
                    ax[axInd].get_ylim()[1],
                    transpTab[['emitt_x', 'emitt_y', 'emitt_4d']].stack().max()])])
            ax[axInd].set_ylim(emitLims/1e3)
            ax[axInd].set_xlabel('s [m]')
            ax[axInd].set_ylabel('Emitt. [pimmrad]')
            ax[axInd].legend(['x (2d)', 'y (2d)', 'Trans. 4d'])
            ax[axInd].grid(True)
        elif quantity == 'Sigmas':
            noMarkers = 20
            markEvery = int(len(s) / noMarkers)
            if markEvery < 1:
                markEvery = 1
            p = ax[axInd].plot(s, transpTab[sigmaXName], '--v', markevery=markEvery)
            ax[axInd].plot(
                s, transpTab[sigmaYName], '--^', markevery=markEvery, color=p[0].get_color())
            ax[axInd].set_xlim(sLims)
            sigmaLims = np.array([
                np.min([
                    ax[axInd].get_ylim()[0], transpTab[[sigmaXName, sigmaYName]].stack().min()]),
                np.max([
                    ax[axInd].get_ylim()[1], transpTab[[sigmaXName, sigmaYName]].stack().max()])])
            ax[axInd].set_ylim(sigmaLims)
            ax[axInd].set_ylabel('Sigma [mm]')
            ax[axInd].grid(True)
        elif quantity == 'TwissBetas':
            try:
                betaLims = np.array([
                    np.min([ax[axInd].get_ylim()[0], transpTab['beta_x'].min()]),
                    np.max([ax[axInd].get_ylim()[1], transpTab['beta_y'].max()])
                ])
                p = ax[axInd].plot(s, transpTab['beta_x'], '--v', markevery=markEvery)
                ax[axInd].plot(
                    s, transpTab['beta_y'], '--^', markevery=markEvery, color=p[0].get_color())
                ax[axInd].set_xlim(sLims)
                ax[axInd].set_ylim(betaLims)
                ax[axInd].set_ylabel('Betatron function [m]')
                ax[axInd].grid(True)
            except KeyError:
                pass
        elif quantity is None:
            ax[axInd].set_visible(False)
        ax[-1].set_xlabel('s [m]')
