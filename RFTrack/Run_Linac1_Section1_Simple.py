import RF_Track as rft
import BeamDynamics as bd
import SimulationData as sd
import numpy as np
import matplotlib.pyplot as plt
import json


### USER INPUT START ###

# BUNCH_TYPE = 'Gaussian'
#
BUNCH_TYPE = 'FromRFTrack'
REL_PATH = 'SimulationRuns/DistrsFromExternalPartners/PositronsAt200MeV/YongkeDistrsV1'
# FILE_NAME = 'CTSB-N02-F100-E06-S0.5-T5.0_FCTest_Pavel_SolC_CLICTW-OptionA08B7-Bc0.50.dat'
# FILTER_SPECS_SELECTOR = ''
FILE_NAME = 'CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_CLICTW-Ztc200-Ri15-Bc0.50.dat'
FILTER_SPECS_SELECTOR = 'CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_CLICTW-Ztc200-Ri15-Bc0.50_MainBunch'
# FILE_NAME = 'CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_PSISW-Ztc200-Ri15-Bc1.50.dat'
# FILTER_SPECS_SELECTOR = 'CTSB-N02-F100-E06-S0.5-T5.0_HTSTest_JNov04_SolC_PSISW-Ztc200-Ri15-Bc1.50_MainBunch'
RFTRACK_FORMAT = 'rftrack_xp_t'

R_APERTURE = 20e-3   # [mm]

SOLENOID_TYPE = 'Analytical'
BZ_SOLENOID = 0.5   # [T]
L_SOLENOID = 0.65   # [m]
R_IN_SOLENOID = 0.100   # [m]
L_SEPARATION = 0.05   # [m]
#
# SOLENOID_TYPE = 'SimulatedFieldmap'
# FILE_PATH = '/afs/psi.ch/project/Pcubed/SimulationRuns/Fieldmaps/custom_solenoid_n_9_cavity.txt'
# BZ_CORR_FACTOR = 0.955   # For L_SOLENOID = 0.65
# BZ_SOLENOID = 0.5   # [T]
# L_SOLENOID = 0.65   # [m]
# L_SEPARATION = 0.05   # [m]

N_SOLENOIDS = 124

EZ_ACC_DRIFT = 15e6   # [V/m]

### USER INPUT END ###


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
    # TODO: Following value to be revised
    Q_DRIVE_BEAM = 1.e-9   # [C]
    N_MACROPARTICLES_DRIVE_BEAM = 10000
    PARTICLE_MASS = rft.electronmass   # [MeV/c/c]
    PARTICLE_CHARGE = +1   # [e], +1 = positrons, -1 = electrons
    # TODO: Remove next line when code is working
    # M0Original = bd.load_rftrack_yongke_1(DISTR_PATH).to_numpy()
    beamIn, _ = bd.convert_rftrack_to_standard_df(
        sourceFilePath=DISTR_PATH, sourceFormat='rftrackYongke1', rftrackDfFormat=RFTRACK_FORMAT,
        z=np.nan, pdgId=-11
    )
    # TODO: Ask exact value to Yongke
    beamIn['z'] = beamIn.loc[0, 't'] * bd.C / 1e6   # [mm]
    beamIn['Q'] = Q_DRIVE_BEAM / N_MACROPARTICLES_DRIVE_BEAM
    with open(sd.build_data_path(REL_PATH, 'filterSpecs.json'), 'r') as filterSpecsFile:
        filterSpecsList = json.load(filterSpecsFile)
    filterSpecs = filterSpecsList[FILTER_SPECS_SELECTOR]['filterSpecs']
    beamIn = bd.filter_distr(beamIn, filterSpecs)
    M0 = bd.convert_standard_df_to_rftrack(standardDf=beamIn, rftrackDfFormat=RFTRACK_FORMAT)[0].to_numpy()
    BUNCH_POPULATION = M0.shape[0]
    B06d = rft.Bunch6d(PARTICLE_MASS, BUNCH_POPULATION, PARTICLE_CHARGE, M0)
    B0 = rft.Bunch6dT(B06d)
print(B0.t)   # Bunch6dT
# print(B0.S)   # Bunch6d

# B = load('../files/beam_S1200.dat.gz').M1;
# B(:,5) -= mean(B(:,5));
t0 = np.mean(B0.get_phase_space('%t0'))
# TODO: Following line not working
# P0 = rft. Bunch6dT([0, 0, 0, 0, 0, 200., rft.electronmass, -1, 0, t0])
# B0 = Bunch6dT(RF_Track.electronmass, 300 * RF_Track.pC, -1, B); 
if SOLENOID_TYPE == 'Analytical':
    zAxis, solBzOnAxis = bd.generate_solenoid_fieldmap(L_SOLENOID, BZ_SOLENOID, R_IN_SOLENOID)
    solFieldmapStep = zAxis[1] - zAxis[0]
elif SOLENOID_TYPE == 'SimulatedFieldmap':
    solMatrix = np.loadtxt(FILE_PATH)
    solFieldmapStep = solMatrix[1,0] - solMatrix[0,0]
    solBzOnAxis = BZ_SOLENOID * BZ_CORR_FACTOR * solMatrix[:,1]
solenoid = rft.Static_Magnetic_FieldMap_1d(solBzOnAxis, solFieldmapStep)
L_SPACING = L_SOLENOID + L_SEPARATION

# Drift
L_DRIFT = N_SOLENOIDS * L_SPACING    # [m]
drift = rft.Drift(L_DRIFT)
drift.set_static_Efield(0, 0, EZ_ACC_DRIFT)

# Tracking volume
vol = rft.Volume()
vol.set_aperture(R_APERTURE)
for solInd in range(0,N_SOLENOIDS):
    vol.add(solenoid, 0, 0, solInd*L_SPACING, 'center')   # element, X, Y, Z in [m]
vol.add(drift, 0, 0, 0)
# Limit tracking
# vol.set_s0(0)   # [m]
# vol.set_s1(10.) # [m]

# Tracking options
TO = rft.TrackingOptions()
TO.dt_mm = 1.   # [mm/c]
TO.odeint_algorithm = 'rk2'   # 'rkf45', 'leapfrog', ...
TO.tt_dt_mm = 10.   # [mm/c], track the emittance every tt_dt_mm steps
# Watch points
# TO.wp_dt_mm = 1e3   # [mm/c], track the emittance every tt_dt_mm steps
# TO.t_max_mm = 1000   # [mm/c]
TO.verbosity = 1   # 0 (default), 1 or 2

# Tracking
# TODO: Uncomment when definition of P0 is working
# finalMomentum = vol.autophase(P0)
B1 = vol.track(B0, TO)
M1 = B1.get_phase_space()   # x [mm], Px [MeV/c], y [mm], Py [MeV/c], S [mm], Pz [MeV/c]
# M1 = B1.get_phase_space('%X %xp %Y %yp %S %P')   # x [mm] Px [MeV/c]   y Py S Pz
# TODO: What was this?
# B1S = vol.get_bunch_at_s1()

# Compute capture efficiency
Mlost = B1.get_lost_particles()   # Columns 1-6 like M1, t [mm/c] at which particle was lost, m [kg], Q [?] of particle type, Q of macro-particle [?]
Mlost = Mlost[Mlost[:, 4].argsort()]
sCapture = Mlost[:, 4]
captureEff = (M0.shape[0] - np.arange(1, Mlost.shape[0]+1, 1)) / M0.shape[0]

TT = vol.get_transport_table('%mean_S %emitt_x %emitt_y %emitt_4d %sigma_X %sigma_Y %mean_E')
beamOut, _ = bd.convert_rftrack_to_standard_df(
    rftrackDf=M1, rftrackDfFormat='rftrack_Px_S', t=B1.t/bd.C*1e6,
    pdgId=-11, Qbunch=Q_DRIVE_BEAM*M1.shape[0]/N_MACROPARTICLES_DRIVE_BEAM
)

# Prepare Bz for plotting
zAxis = np.linspace(vol.get_s0(), vol.get_s1(), 1000)   # [m]
Bz = []
for z in zAxis:
    [E, B] = vol.get_field(0, 0, z*1e3, 0)   # x,y,z,t (mm, mm/c)
    Bz.append(B[2])
Bz = np.array(Bz)

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

sLims = TT[[0, -1], 0] / 1e3
Elims = [0, TT[-1,6]]
fig2, ax2 = plt.subplots(5, 1)
ax2[0].plot(zAxis, Bz)
ax2[0].set_xlim(sLims)
ax2[0].set_xlabel('s [m]')
ax2[0].set_ylabel('Bz [T]')
ax2[0].grid()
ax2[1].plot(TT[:,0]/1e3, TT[:,6])
ax2[1].set_xlim(sLims)
ax2[1].set_ylim(Elims)
ax2[1].set_xlabel('s [m]')
ax2[2].set_ylabel('Beam energy [MeV]')
ax2[2].plot(sCapture/1e3, captureEff)
ax2[2].set_xlim(sLims)
ax2[2].set_ylim([0, 1])
ax2[2].set_xlabel('s [m]')
ax2[2].set_ylabel('Capture efficiency')
ax2[2].grid()
ax2[3].plot(TT[:,0]/1e3, TT[:,1])
ax2[3].plot(TT[:,0]/1e3, TT[:,2])
ax2[3].plot(TT[:,0]/1e3, TT[:,3])
ax2[3].set_xlim(sLims)
ax2[3].set_xlabel('s [m]')
ax2[3].set_ylabel('emitt [pi mm mrad]')
ax2[3].legend(['emitt x', 'emitt y', 'emitt 4d'])
ax2[3].grid()
ax2[4].plot(TT[:,0]/1e3, TT[:,4])
ax2[4].plot(TT[:,0]/1e3, TT[:,5])
ax2[4].set_xlim(sLims)
ax2[4].set_xlabel('s [m]')
ax2[4].set_ylabel('sigma [mm]')
ax2[4].legend(['sigma x', 'sigma y'])
ax2[4].grid()

# # plt.ion()
plt.show(block=False)
input("Press Enter to continue...")
# plt.close('all')
