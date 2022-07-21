import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import BeamDynamics as bd
import RF_Track as rft


DEFAULT_COLOR_CYCLE = plt.rcParams["axes.prop_cycle"].by_key()['color']

# TODO: Make path relative
RF_CLIC_FIELDMAP_BASEPATH = \
    '/home/tia/Repos/GIT_PSIPositronProduction/' + \
    'RFTrack/YongkeTool_V1/field/field_map_CLIC_Lband'
RF_CLIC_FREQ = 1.9986163867e+09   # [Hz]
RF_CLIC_GRADIENT = 11.23e6   # [V/m]
RF_CLIC_CELLS_PER_PERIOD = 3.
RF_CLIC_L_CELL = bd.C / RF_CLIC_FREQ / RF_CLIC_CELLS_PER_PERIOD  # [m]


def load_fieldmap_rf_clic():
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
        rfFields[complComp] = rfFields[complComp].reshape(matDimensions) \
            .transpose()   # .astype(np.complex128)
    return rfFields, meshes


def rf_clic_single_period(rfFieldmapDim):
    # Fields of 1 RF period (3 cells)
    # TODO: Data location and path construction
    rfFields, rfMesh = load_fieldmap_rf_clic()
    if rfFieldmapDim == '1D':
        rf = rft.RF_FieldMap_1d(
            rfFields['Ez'][0, 0, :],
            rfMesh['za'][1] - rfMesh['za'][0],
            rfMesh['za'][-1],
            RF_CLIC_FREQ, +1
        )
    elif rfFieldmapDim == '3D':
        rf = rft.RF_FieldMap(
            rfFields['Ex'], rfFields['Ey'], rfFields['Ez'],
            rfFields['Bx'], rfFields['By'], rfFields['Bz'],
            rfMesh['ra'][0], rfMesh['ta'][0],
            rfMesh['ra'][1] - rfMesh['ra'][0],
            rfMesh['ta'][1] - rfMesh['ta'][0],
            rfMesh['za'][1] - rfMesh['za'][0],
            rfMesh['za'][-1],
            RF_CLIC_FREQ, +1
        )
        rf.set_cylindrical(True)
    else:
        raise ValueError('Invalid rfFieldmapDim={:s}.'.format(rfFieldmapDim))
    return rf


def save_plot_transport(ax, vol, beam0, beam1, outRelPath):
    # Prepare Ez and Bz for plotting
    zAxis = np.linspace(vol.get_s0(), vol.get_s1(), 1000)   # [m]
    Ez = []
    Bz = []
    for z in zAxis:
        E, B = vol.get_field(0, 0, z*1e3, 0)   # x,y,z,t (mm, mm/c)
        Ez.append(E[2])
        Bz.append(B[2])
    Ez = np.array(Ez)
    Bz = np.array(Bz)
    emFields = pd.DataFrame(np.row_stack([zAxis, Ez, Bz]).T, columns=['z', 'Ez', 'Bz'])
    emFields.to_csv(os.path.join(outRelPath, 'EMFields.dat'), index=None)
    # Get transport table
    getTransportTableStr = '%mean_S %emitt_x %emitt_y %emitt_4d %sigma_X %sigma_Y %mean_E'
    TT = vol.get_transport_table(getTransportTableStr)
    transportTable = pd.DataFrame(TT, columns=getTransportTableStr.replace('%', '').split())
    transportTable.to_csv(os.path.join(outRelPath, 'TransportTable.dat'), index=None)
    # Compute capture efficiency
    M0 = beam0.get_phase_space()
    Mlost = beam1.get_lost_particles()
    # Columns of Mlost like columns 1-6 of M0, in addition:
    # t [mm/c] at which particle was lost, m [kg],
    # Q [?] of particle type, Q of macro-particle [?]
    Mlost = Mlost[Mlost[:, 4].argsort()]
    sCapture = Mlost[:, 4]
    captureEff = (
        M0.shape[0] - np.arange(1, Mlost.shape[0]+1, 1)
    ) / M0.shape[0]
    captureEfficiency = pd.DataFrame(
        np.row_stack([sCapture, captureEff]).T, columns=['s', 'CaptureEfficiency']
    )
    captureEfficiency.to_csv(os.path.join(outRelPath, 'CaptureEfficiency.dat'), index=None)
    plot_transport(ax, emFields, transportTable, captureEfficiency)


def load_plot_transport(ax, simRelPath, sShiftEMFields=0, tShift=0):
    # TODO: Global variable, use also in save_plot_transport()
    transportFileNames = ['EMFields.dat', 'TransportTable.dat', 'CaptureEfficiency.dat']
    transportDfs = []
    for fileName in transportFileNames:
        transportDfs.append(
            pd.read_csv(os.path.join(simRelPath, fileName))
        )
    plot_transport(ax, *transportDfs, sShiftEMFields=sShiftEMFields, tShift=tShift)


# TODO: Move following function to module BeamDynamics?
def plot_transport(ax, emFields, transpTab, captureEff, sShiftEMFields=0, tShift=0):
    try:
        s = transpTab['mean_S'] / 1e3  # [m]
        sigmaXName = 'sigma_X'
        sigmaYName = 'sigma_Y'
    except KeyError:
        s = (transpTab['mean_t'] + tShift) / 1e3  # [m/c]
        sigmaXName = 'sigma_x'
        sigmaYName = 'sigma_y'
    # TODO: Recognize when axis is empty (xlim = [0, 1])
    sLims = np.array([np.min([ax[0].get_xlim()[0], s.min()]), np.max([ax[0].get_xlim()[1], s.max()])])
    BzLims = np.array([0, np.max([ax[0].get_ylim()[1], emFields['Bz'].max()])])
    EzLims = np.array([
        np.min([ax[1].get_ylim()[0], emFields['Ez'].min()]),
        np.max([ax[1].get_ylim()[1], emFields['Ez'].max()])
    ])
    Elims = np.array([
        np.min([ax[2].get_ylim()[0], transpTab['mean_E'].min()]),
        np.max([ax[2].get_ylim()[1], transpTab['mean_E'].max()])
    ])
    emitLims = np.array([
        np.min([ax[4].get_ylim()[0], transpTab[['emitt_x', 'emitt_y', 'emitt_4d']].stack().min()]),
        np.max([ax[4].get_ylim()[1], transpTab[['emitt_x', 'emitt_y', 'emitt_4d']].stack().max()])
    ])
    sigmaLims = np.array([
        np.min([ax[5].get_ylim()[0], transpTab[[sigmaXName, sigmaYName]].stack().min()]),
        np.max([ax[5].get_ylim()[1], transpTab[[sigmaXName, sigmaYName]].stack().max()])
    ])
    ax[0].plot(emFields['z']+sShiftEMFields, emFields['Bz'])  # , color=DEFAULT_COLOR_CYCLE[0]
    ax[0].set_xlim(sLims)
    ax[0].set_ylim(BzLims)
    ax[0].set_xlabel('s [m]')
    ax[0].set_ylabel('Bz [T]')  # , color=DEFAULT_COLOR_CYCLE[0]
    # ax[0].grid()
    ax[1].plot(emFields['z']+sShiftEMFields, emFields['Ez']/1e6, '--')  # , color=DEFAULT_COLOR_CYCLE[0]
    ax[1].set_ylim(EzLims/1e6)
    ax[1].set_ylabel('Ez [MV/m]')  # , color=DEFAULT_COLOR_CYCLE[0]
    ax[2].plot(s, transpTab['mean_E'])
    ax[2].set_xlim(sLims)
    ax[2].set_ylim(Elims)
    ax[2].set_xlabel('s [m]')
    ax[2].set_ylabel('E [MeV]')
    # ax[2].grid()
    ax[3].plot(captureEff['s']/1e3+sShiftEMFields, captureEff['CaptureEfficiency'], '--')
    ax[3].set_ylim([0, 1])
    ax[3].set_ylabel('Capture eff.')
    noMarkers = 20
    markEvery = int(len(s) / noMarkers)
    p = ax[4].plot(s, transpTab['emitt_x']/1e3, '-v', markevery=markEvery)
    ax[4].plot(s, transpTab['emitt_y']/1e3, '-^', markevery=markEvery, color=p[0].get_color())
    ax[4].plot(s, transpTab['emitt_4d']/1e3, '-', markevery=markEvery, color=p[0].get_color())
    ax[4].set_xlim(sLims)
    ax[4].set_ylim(emitLims/1e3)
    ax[4].set_xlabel('s [m]')
    ax[4].set_ylabel('Emitt. [pimmrad]')
    ax[4].legend(['x (2d)', 'y (2d)', 'Trans. 4d'])
    # ax[4].grid()
    p = ax[5].plot(s, transpTab[sigmaXName], '--v', markevery=markEvery)
    ax[5].plot(s, transpTab[sigmaYName], '--^', markevery=markEvery, color=p[0].get_color())
    ax[5].set_ylim(sigmaLims)
    ax[5].set_ylabel('Sigma [mm]')
