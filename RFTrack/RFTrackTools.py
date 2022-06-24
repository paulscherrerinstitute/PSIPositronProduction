import numpy as np
import BeamDynamics as bd
import RF_Track as rft


# TODO: Make path relative
RF_CLIC_FIELDMAP_BASEPATH = \
    '/home/tia/Repos/GIT_PSIPositronProduction/' + \
    'RFTrack/YongkeTool_V1/field/field_map_CLIC_Lband'
RF_CLIC_FREQ = 1.9986163867e+09   # [Hz]
RF_CLIC_GRADIENT = 11.23e6   # [V/m]
RF_CLIC_CELLS_PER_PERIOD = 3.
RF_CLIC_L_CELL = bd.C / RF_CLIC_FREQ / RF_CLIC_CELLS_PER_PERIOD


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


def plot_transport(ax, vol, beam0, beam1):
    # Prepare Bz for plotting
    zAxis = np.linspace(vol.get_s0(), vol.get_s1(), 1000)   # [m]
    Bz = []
    for z in zAxis:
        E, B = vol.get_field(0, 0, z*1e3, 0)   # x,y,z,t (mm, mm/c)
        Bz.append(B[2])
    Bz = np.array(Bz)
    # Get transport table
    TT = vol.get_transport_table(
        '%mean_S %emitt_x %emitt_y %emitt_4d %sigma_X %sigma_Y %mean_E'
    )
    sLims = np.array([np.min(TT[:, 0]), np.max(TT[:, 0])]) / 1e3
    Elims = [0, np.max(TT[:, 6])]
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
    # Plot
    ax[0].plot(zAxis, Bz)
    ax[0].set_xlim(sLims)
    ax[0].set_xlabel('s [m]')
    ax[0].set_ylabel('Bz [T]')
    ax[0].grid()
    ax[1].plot(TT[:, 0]/1e3, TT[:, 6])
    ax[1].set_xlim(sLims)
    ax[1].set_ylim(Elims)
    ax[1].set_xlabel('s [m]')
    ax[1].set_ylabel('Beam energy [MeV]')
    ax[1].grid()
    ax[2].plot(sCapture/1e3, captureEff)
    ax[2].set_xlim(sLims)
    ax[2].set_ylim([0, 1])
    ax[2].set_xlabel('s [m]')
    ax[2].set_ylabel('Capture efficiency')
    ax[2].grid()
    ax[3].plot(TT[:, 0]/1e3, TT[:, 1])
    ax[3].plot(TT[:, 0]/1e3, TT[:, 2])
    ax[3].plot(TT[:, 0]/1e3, TT[:, 3])
    ax[3].set_xlim(sLims)
    ax[3].set_xlabel('s [m]')
    ax[3].set_ylabel('emitt [pi mm mrad]')
    ax[3].legend(['emitt x', 'emitt y', 'emitt 4d'])
    ax[3].grid()
    ax[4].plot(TT[:, 0]/1e3, TT[:, 4])
    ax[4].plot(TT[:, 0]/1e3, TT[:, 5])
    ax[4].set_xlim(sLims)
    ax[4].set_xlabel('s [m]')
    ax[4].set_ylabel('sigma [mm]')
    ax[4].legend(['sigma x', 'sigma y'])
    ax[4].grid()
