import RF_Track as rft
import RFTrackTools as rfttools
import BeamDynamics as bd
import SimulationData as sd
import numpy as np
import matplotlib.pyplot as plt
import json


# USER INPUT START ###

# BUNCH_TYPE = 'Gaussian'
#
BUNCH_TYPE = 'FromRFTrack'
REL_PATH = 'SimulationRuns/DistrsFromExternalPartners/PositronsAt200MeV/YongkeDistrsV1'
# FILTER_SPECS_SELECTOR = None
# FILE_NAME = 'CTSB-N02-F100-E06-S0.5-T5.0_FCTest_Pavel_SolC_CLICTW-OptionA08B7-Bc0.50.dat'
# FILTER_SPECS_SELECTOR = ''
FILE_NAME = 'CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_CLICTW-Ztc200-Ri15-Bc0.50.dat'
FILTER_SPECS_SELECTOR = \
    'CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_CLICTW-Ztc200-Ri15-Bc0.50_MainBunch'
# FILE_NAME = 'CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_PSISW-Ztc200-Ri15-Bc1.50.dat'
# FILTER_SPECS_SELECTOR = \
#     'CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_PSISW-Ztc200-Ri15-Bc1.50_MainBunch'
RFTRACK_FORMAT = 'rftrack_xp_t'

R_APERTURE = 20e-3   # [m]

# SOLENOID_TYPE = 'HomogeneousChannel'
# BZ_SOLENOID = 0.5   # [T]
# L_SOLENOID = 0.65   # [m]
# SEPARATION_SOLENOID = 0.05   # [m]
#
SOLENOID_TYPE = 'Analytical'
BZ_SOLENOID = 0.5   # [T]
R_IN_SOLENOID = 0.100   # [m]
SEPARATION_SOLENOID = 0.3   # [m]
# TODO: Give better names to the two following variables
FILE_PATH = None
BZ_CORR_FACTOR = None
#
# SOLENOID_TYPE = 'SimulatedFieldmap'
# FILE_PATH = '/afs/psi.ch/project/Pcubed/SimulationRuns/Fieldmaps/custom_solenoid_n_9_cavity.txt'
# BZ_CORR_FACTOR = 0.955   # For L_SOLENOID = 0.65
# BZ_SOLENOID = 0.5   # [T]
# L_SOLENOID = 0.65   # [m]
# SEPARATION_SOLENOID = 0.05   # [m]

# E_GAIN_LINAC = 0.06e9   # [eV]

# RF_TYPE = 'HomogeneousAccel'
# RF_GRADIENT = 15e6   # [V/m]
#
AUTOPHASING = False
RF_TYPE = 'field_map_CLIC_Lband'
RF_FIELDMAP_DIM = '3D'
RF_N_CELLS = 30
RF_N_PERIODS_PER_STRUCTURE = np.round(RF_N_CELLS/3)
RF_L_STRUCTURE = RF_N_PERIODS_PER_STRUCTURE * bd.C / rfttools.RF_CLIC_FREQ   # [m]
RF_N_STRUCTURES = 2  # 58
RF_GRADIENT = 15.0E6   # [V/m]
RF_T0 = -70./360. / rfttools.RF_CLIC_FREQ * bd.C * 1e3   # [mm/c]
RF_PHASE = 0   # [deg]
RF_SEPARATION = 0.2   # [m]

RF_HARMONIC_ON = False
RF_HARMONIC_FREQ = 3 * rfttools.RF_CLIC_FREQ
RF_HARMONIC_PHASE_ADV = 2./3.*np.pi   # [rad]
RF_HARMONIC_N_CELLS = 120
RF_HARMONIC_GRADIENT = 30e6   # [V/m]
RF_HARMONIC_PHASE = 30.0   # [deg]

L_SOLENOID = RF_L_STRUCTURE + RF_SEPARATION - SEPARATION_SOLENOID
SPACING_SOLENOID = L_SOLENOID + SEPARATION_SOLENOID
N_SOLENOIDS = RF_N_STRUCTURES + 1
if RF_HARMONIC_ON:
    N_SOLENOIDS += 1
# N_SOLENOIDS = int(E_GAIN_LINAC/RF_GRADIENT/SPACING_SOLENOID)

# USER INPUT END ###


if BUNCH_TYPE == 'Gaussian':
    Q = 4.   # [nC]
    N_MACROPARTICLES = 10000   # Macro-particles
    SIGMA_Z = 0.1   # [mm]
    P_SPREAD = 100.   # [permil], momentum spread
    P_REF = 200.   # [MeV/c]
    Twiss = rft.Bunch6dT_twiss()
    Twiss.emitt_x = 10e3   # [pi mm mrad]
    Twiss.emitt_y = 10e3   # [pi mm mrad]
    Twiss.beta_x = 1.   # [m]
    Twiss.beta_y = 1.   # [m]
    Twiss.alpha_x = 0.01
    Twiss.alpha_y = 0.01
    Twiss.sigma_z = SIGMA_Z   # [mm/c]
    Twiss.sigma_d = P_SPREAD   # [permil]
    B0 = rft.Bunch6dT(
        rft.electronmass, Q*rft.nC, +1, P_REF, Twiss, N_MACROPARTICLES
    )
    M0 = B0.get_phase_space()   # x [mm], Px [MeV/c], y [mm], Py [MeV/c], S [mm], Pz [MeV/c]
    beamIn, _ = bd.convert_rftrack_to_standard_df(
        rftrackDf=M0, rftrackDfFormat='rftrack_Px_S', t=np.nan, pdgId=-11, Qbunch=Q*1e-9
    )
elif BUNCH_TYPE == 'FromRFTrack':
    DISTR_PATH = sd.build_data_path(REL_PATH, FILE_NAME)
    # TODO: Following value to be reviewed
    Q_DRIVE_BEAM = 1.e-9   # [C]
    N_MACROPARTICLES_DRIVE_BEAM = 10000
    PARTICLE_MASS = rft.electronmass   # [MeV/c/c]
    PARTICLE_CHARGE = +1   # [e], +1 = positrons, -1 = electrons
    beamIn, _ = bd.convert_rftrack_to_standard_df(
        sourceFilePath=DISTR_PATH, sourceFormat='rftrackYongke1', rftrackDfFormat=RFTRACK_FORMAT,
        z=np.nan, pdgId=-11
    )
    # TODO: Ask exact value to Yongke
    beamIn['z'] = beamIn.loc[0, 't'] * bd.C / 1e6   # [mm]
    beamIn['Q'] = Q_DRIVE_BEAM / N_MACROPARTICLES_DRIVE_BEAM
    if FILTER_SPECS_SELECTOR is not None:
        with open(sd.build_data_path(REL_PATH, 'filterSpecs.json'), 'r') as filterSpecsFile:
            filterSpecsList = json.load(filterSpecsFile)
        filterSpecs = filterSpecsList[FILTER_SPECS_SELECTOR]['filterSpecs']
        refParticle = filterSpecsList[FILTER_SPECS_SELECTOR]['RefParticle1']
        beamIn = bd.filter_distr(beamIn, filterSpecs)
    M0 = bd.convert_standard_df_to_rftrack(
        standardDf=beamIn, rftrackDfFormat=RFTRACK_FORMAT
    )[0].to_numpy()
    BUNCH_POPULATION = M0.shape[0]
    B06d = rft.Bunch6d(PARTICLE_MASS, BUNCH_POPULATION, PARTICLE_CHARGE, M0)
    B0 = rft.Bunch6dT(B06d)

# TODO: Document following 2 lines
# print(B0.t)   # Bunch6dT
# print(B0.S)   # Bunch6d

# B = load('../files/beam_S1200.dat.gz').M1;
# B(:,5) -= mean(B(:,5));
T0_REF = 62.867e-9 * bd.C * 1e3   # [mm/c]
E0_REF = 184.
# TODO: other option could be: np.mean(B0.get_phase_space('%t0'))
# TODO: Check that particle charge is correct in the following line
P0 = rft. Bunch6dT(np.array([
    refParticle['x'],
    refParticle['xp'],
    refParticle['y'],
    refParticle['yp'],
    0,   # z
    refParticle['pz'],
    rft.electronmass,
    PARTICLE_CHARGE,
    0,
    refParticle['t'] * bd.C / 1e6   # [mm/c]
]))
# B0 = Bunch6dT(RF_Track.electronmass, 300 * RF_Track.pC, -1, B);


if SOLENOID_TYPE == 'Analytical':
    zAxis, solBzOnAxis = bd.generate_solenoid_fieldmap(L_SOLENOID, BZ_SOLENOID, R_IN_SOLENOID)
    solFieldmapStep = zAxis[1] - zAxis[0]
elif SOLENOID_TYPE == 'SimulatedFieldmap':
    solMatrix = np.loadtxt(FILE_PATH)
    solFieldmapStep = solMatrix[1, 0] - solMatrix[0, 0]
    solBzOnAxis = BZ_SOLENOID * BZ_CORR_FACTOR * solMatrix[:, 1]


if RF_TYPE == 'HomogeneousAccel':
    L_DRIFT = N_SOLENOIDS * SPACING_SOLENOID    # [m]
    drift = rft.Drift(L_DRIFT)
    drift.set_static_Efield(0, 0, RF_GRADIENT)
    if SOLENOID_TYPE == 'HomogeneousChannel':
        drift.set_static_Bfield(0, 0, BZ_SOLENOID)
else:
    rf = rfttools.rf_clic_single_period(RF_FIELDMAP_DIM)
    rfGap = rft.Drift(RF_SEPARATION)
    if SOLENOID_TYPE == 'HomogeneousChannel':
        rf.set_static_Bfield(0, 0, BZ_SOLENOID)
        rfGap.set_static_Bfield(0, 0, BZ_SOLENOID)
    lat = rft.Lattice()
    # Initial gap
    INITIAL_RF_GAP_L = (L_SOLENOID + SEPARATION_SOLENOID + RF_SEPARATION) / 2.
    initialRfGap = rft.Drift(INITIAL_RF_GAP_L)
    tRf = RF_T0 - INITIAL_RF_GAP_L*1e3
    # TODO: Check double drift space at the beginning
    lat.append(initialRfGap)
    for structInd in np.arange(RF_N_STRUCTURES):
        lat.append(rfGap)
        rf.set_phid(RF_PHASE)
        rf.set_P_actual((RF_GRADIENT/rfttools.RF_CLIC_GRADIENT)**2.)
        rf.set_t0(tRf)
        for rfPeriodInd in np.arange(RF_N_PERIODS_PER_STRUCTURE):
            lat.append(rf)
        tRf += RF_SEPARATION * 1e3   # [mm/c]

if RF_HARMONIC_ON:
    # TODO: Document this
    # By N_CELLS:
    # - = TW starts at the beginning of the cell
    # + = TW starts in the middle of the cell (useful for couplers)
    RF_HARMONIC_MAIN_MODE_INDEX = 0
    rfHarmonic = rft.TW_Structure(
        np.array([RF_HARMONIC_GRADIENT]), RF_HARMONIC_MAIN_MODE_INDEX,
        RF_HARMONIC_FREQ, RF_HARMONIC_PHASE_ADV, -1.*RF_HARMONIC_N_CELLS
    )
    rfHarmonic.set_t0(tRf)
    rfHarmonic.set_phid(RF_HARMONIC_PHASE)
    lat.append(rfGap)
    lat.append(rfHarmonic)


# Tracking volume
vol = rft.Volume()
# TODO: Is this a circular aperture?
vol.set_aperture(R_APERTURE)
if RF_TYPE == 'HomogeneousAccel':
    vol.add(drift, 0, 0, 0)
else:
    vol.add(lat, 0, 0, 0)
if not SOLENOID_TYPE == 'HomogeneousChannel':
    solenoid = rft.Static_Magnetic_FieldMap_1d(solBzOnAxis, solFieldmapStep)
    for solInd in range(0, N_SOLENOIDS):
        vol.add(solenoid, 0, 0, solInd*SPACING_SOLENOID, 'center')   # element, X, Y, Z in [m]

print('Setup complete, going to track...')

# Limit tracking
# vol.set_s0(0)   # [m]
# vol.set_s1(10.) # [m]

# Tracking options
TO = rft.TrackingOptions()
TO.dt_mm = 1.   # [mm/c]
TO.odeint_algorithm = 'rk2'   # 'rkf45', 'leapfrog', ...
TO.tt_dt_mm = 10.   # [mm/c], track the emittance every tt_dt_mm (time)
# Watch points
# TO.wp_dt_mm = 1e3   # [mm/c], save the distribution on disk every wp_dt_mm (time)
# TO.t_max_mm = 1000   # [mm/c]
TO.verbosity = 1   # 0 (default), 1 or 2

# Tracking
if AUTOPHASING:
    finalMomentum = vol.autophase(P0)
    print('Final momentum after autophasing: {:.1f} MeV/c'.format(finalMomentum))
B1 = vol.track(B0, TO)
M1 = B1.get_phase_space()   # x [mm], Px [MeV/c], y [mm], Py [MeV/c], S [mm], Pz [MeV/c]
# M1 = B1.get_phase_space('%X %xp %Y %yp %S %P')   # x [mm] Px [MeV/c]   y Py S Pz
# TODO: What was this?
# B1S = vol.get_bunch_at_s1()

beamOut, _ = bd.convert_rftrack_to_standard_df(
    rftrackDf=M1, rftrackDfFormat='rftrack_Px_S', t=B1.t/bd.C*1e6,
    pdgId=-11, Qbunch=Q_DRIVE_BEAM*M1.shape[0]/N_MACROPARTICLES_DRIVE_BEAM,
    outFwfPath='./DistrOut_Linac1_Section1_Simple'
)

plotDefs = [
    {
        'varName1': 'x', 'varName2': 'y',
        'opacityHist': 0.6,
    },
    {
        'varName1': 'x', 'varName2': 'xp',
        'opacityHist': 0.6,
    },
    {
        'varName1': 'y', 'varName2': 'yp',
        'opacityHist': 0.6,
    },
    {
        'varName1': 'z', 'varName2': 'pz',
        'var1Center': True,
        'opacityHist': 0.6,
    },
    {
        'varName1': 'z', 'varName2': 'Ekin',
        'var1Center': True,
        'opacityHist': 0.6,
    },
]
ax = bd.plot_distr([beamIn, beamOut], plotDefs)

fig2, ax2 = plt.subplots(5, 1)
rfttools.plot_transport(ax2[1:], vol, B0, B1)

# # plt.ion()
plt.show(block=False)
input("Press Enter to continue...")
# plt.close('all')
