import RF_Track as rft
import RFTrackTools as rfttools
import BeamDynamics as bd
import SimulationData as sd
import numpy as np
import matplotlib.pyplot as plt
import json


BUNCH_TYPE = 'FromGEANT4'
BUNCH_FILEPATH = 'RFTrack/YongkeTool_V2/amd_input/E6GeV_SpotSize0.5mm_Target5X0.dat'
RFTRACK_FORMAT = 'rftrack_xp_t'
BUNCH_PDGID = -11
PARTICLE_CHARGE = +1   # [e], +1 = positrons, -1 = electrons
PARTICLE_MASS = rft.electronmass   # [MeV/c/c]
BUNCH_Z = 0.
BUNCH_DOWNSAMPLING = 1

TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = +41e-3   # [m]
#
# RF_TYPE = 'field_map_PSI_Sband'
# TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = +20e-3   # [m]
# SOL_HOMOG_BZ = 1.5   # [T]

# TODO
# target.Ne = 1e4;	% no. of e- simulated
# target.Np = 4.37e10;	% no. of e+ required

VOL_R_APERTURE = 1.  # [m]

AMD_FIELDMAP_2D = 'RFTrack/YongkeTool_V2/field/field_map_HTS_5coils_Apr2022.dat'
AMD_FIELDMAP_1D = 'RFTrack/YongkeTool_V2/field/field_map_HTS_5coils_Apr2022_1D.dat'
USE_AMD_FIELDMAP_3D = True
AMD_R_APERTURE = 20e-3   # [m]
# TODO: This is an old value
AMD_L_HALF_MECHANICAL = 96.5e-3   # [m]

TRACK_AFTER_AMD = False

RF_FIELDMAP_DIM = '3D'
RF_N_STRUCTURES = 11  # 11
RF_N_CELLS = 30
RF_N_PERIODS_PER_STRUCTURE = np.round(RF_N_CELLS/rfttools.RF_CLIC_CELLS_PER_PERIOD)
RF_L_STRUCTURE = RF_N_PERIODS_PER_STRUCTURE * bd.C/rfttools.RF_CLIC_FREQ  # [m]
RF_L_FLANGE = 0.05461  # [m]
RF_L_MECH_MARGIN = 0.01  # [m]
RF_T0 = 17.5 + 287.2  # [mm/c]
RF_PHASE_1ST = 171.  # [deg]
RF_PHASE_AFTER_1ST = 171.  # [deg]
RF_GRADIENT_1ST = 17.5e6  # [V/m]
RF_GRADIENT_AFTER_1ST = 21.0e6   # [V/m]
INITIAL_L = 0.   # [m]
RF_SEPARATION = 0.2   # [m]
RF_R_APERTURE = 20e-3   # [m]

# SOLENOID_TYPE = 'HomogeneousChannel'
# SOL_HOMOG_BZ = 0.5   # [T]
#
SOLENOID_TYPE = 'Analytical'
SOL_R_IN_COIL = 0.130   # [m]
SOL_R_OUT_COIL = 0.250   # [m]
SOL_L = RF_L_STRUCTURE - RF_L_FLANGE - rfttools.RF_CLIC_L_CELL - RF_L_MECH_MARGIN
SOL_J = 3.54e6   # [A/m2]
SOL_HOMOG_BZ = 0.5   # [T]

FINAL_L = 0.   # [m]
T_ADD_NON_RELATIVISTIC = np.Inf  # 1000.  # [mm/c]
N_RF_STRUCT_1ST_TRACKING = RF_N_STRUCTURES
N_RF_CELLS_LONG_PS_CUT = 1 + 2 * rfttools.RF_CLIC_CELLS_PER_PERIOD
# E.g.: 1 + 2 * ... = keep 3 positron buckets

splitTracking = RF_N_STRUCTURES > N_RF_STRUCT_1ST_TRACKING

OUT_REL_PATH = './Results_AMD/LatestSim/'

# TODO
# rf.R1i  = 20; % mm

# TODO: Refactorize, same code in Run_Linac1_Section1_Simple.py
# A_target = load(BUNCH_FILEPATH)
distrMatNp = np.loadtxt(BUNCH_FILEPATH, skiprows=1)
beamIn, _ = bd.convert_rftrack_to_standard_df(
    rftrackDf=distrMatNp[::BUNCH_DOWNSAMPLING, :], rftrackDfFormat=RFTRACK_FORMAT,
    z=BUNCH_Z, pdgId=BUNCH_PDGID
)
# TODO:
# beamIn['Q'] = Q_DRIVE_BEAM / N_MACROPARTICLES_DRIVE_BEAM
# TODO:
# M = A_target(:,5) < 1000; % mm/c, remove very large time particles
# A_target = A_target(M,:);
# if FILTER_SPECS_SELECTOR is not None:
#     with open(sd.build_data_path(REL_PATH, 'filterSpecs.json'), 'r') \
#         as filterSpecsFile:
#         filterSpecsList = json.load(filterSpecsFile)
#     filterSpecs = filterSpecsList[FILTER_SPECS_SELECTOR]['filterSpecs']
#     refParticle = filterSpecsList[FILTER_SPECS_SELECTOR]['RefParticle1']
#     beamIn = bd.filter_distr(beamIn, filterSpecs)
M0 = bd.convert_standard_df_to_rftrack(
    standardDf=beamIn, rftrackDfFormat=RFTRACK_FORMAT
)[0].to_numpy()
BUNCH_POPULATION = M0.shape[0]
B0_6d = rft.Bunch6d(PARTICLE_MASS, BUNCH_POPULATION, PARTICLE_CHARGE, M0)
B0_6dT = rft.Bunch6dT(B0_6d)

TARGET_EXIT_Z_IN_VOLUME = 0.   # [m]
vol = rft.Volume()
vol.set_aperture(VOL_R_APERTURE, VOL_R_APERTURE, 'circular')

if USE_AMD_FIELDMAP_3D:
    amdFieldmap = bd.load_octave_matrices(AMD_FIELDMAP_2D)
    amdDz = amdFieldmap['Z'][0, 1] - amdFieldmap['Z'][0, 0]  # [mm]
    amdDr = amdFieldmap['R'][1, 0]-amdFieldmap['R'][0, 0]  # [mm]
    amdZ = amdFieldmap['Z'][0, :]  # [mm]
    amdBzOnAxis = amdFieldmap['Bz'][0, :]
else:
    amdFieldmap = bd.load_octave_matrices(AMD_FIELDMAP_1D)
    amdDz = amdFieldmap['Z'][1] - amdFieldmap['Z'][0]  # [mm]
    amdZ = amdFieldmap['Z']  # [mm]
    amdBzOnAxis = amdFieldmap['Bz']

# Get effective field and length used in tracking,
# starts from target exit, stops at constant solenoid field value
indZTargetExit = (np.abs(amdZ - TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD*1e3)).argmin()
# TODO: Make positioning of 1st structure more clear
# (avoid use of SOL_HOMOG_BZ when using analytical fieldmap?)
indZConstBzStart = np.nonzero(amdBzOnAxis > SOL_HOMOG_BZ)[0][-1] + 1
amdFieldLength = amdZ[indZConstBzStart] - amdZ[indZTargetExit]
if SOLENOID_TYPE == 'HomogeneousChannel':
    indsEff = np.arange(indZTargetExit, indZConstBzStart+1)
elif SOLENOID_TYPE == 'Analytical':
    indsEff = np.arange(indZTargetExit, amdBzOnAxis.shape[0])
    amdFieldLengthAnalytical = \
        amdZ[amdBzOnAxis.shape[0]-1] - amdZ[indZTargetExit]
amdZEff = amdZ[indsEff]
amdExitZInVolume = TARGET_EXIT_Z_IN_VOLUME + amdFieldLength*1e-3   # [m]
print('amdFieldLength = {:f} mm.'.format(amdFieldLength))
if USE_AMD_FIELDMAP_3D:
    amdBrEff = amdFieldmap['Br'][:, indsEff]
    amdBzEff = amdFieldmap['Bz'][:, indsEff]
    amd = rft.Static_Magnetic_FieldMap_2d(amdBrEff.T, amdBzEff.T, amdDr*1e-3, amdDz*1e-3)
else:
    amdBzEff = amdFieldmap['Bz'][indsEff]
    amd = rft.Static_Magnetic_FieldMap_1d(amdBzEff.T, amdDz*1e-3)
if SOLENOID_TYPE == 'HomogeneousChannel':
    # TODO: Why is set_length() necessary?
    amd.set_length(amdFieldLength*1e-3)
    # TODO: Parametrize number of steps, i.e. following factor 2
    # TODO: Is the following line necessary in volume?
    amd.set_nsteps(int(amdFieldLength*2.))
elif SOLENOID_TYPE == 'Analytical':
    amd.set_length(amdFieldLengthAnalytical*1e-3)
    # TODO: Is the following line necessary in volume?
    amd.set_nsteps(int(amdFieldLengthAnalytical*2.))
# amd.set_odeint_algorithm('rkf45')
vol.add(amd, 0, 0, TARGET_EXIT_Z_IN_VOLUME, 'entrance')
zFinalInVolume = amdExitZInVolume

AMD_MECH_EXIT_Z_IN_VOLUME = \
    TARGET_EXIT_Z_IN_VOLUME - TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD + AMD_L_HALF_MECHANICAL
amdAperture = rft.Drift(AMD_MECH_EXIT_Z_IN_VOLUME)
amdAperture.set_aperture(AMD_R_APERTURE, AMD_R_APERTURE, 'circular')
vol.add(amdAperture, 0, 0, 0, 'entrance')

if TRACK_AFTER_AMD:
    rf = rfttools.rf_clic_single_period(RF_FIELDMAP_DIM)
    rf.set_aperture(RF_R_APERTURE, RF_R_APERTURE, 'circular')
    lat = rft.Lattice()
    if INITIAL_L > 0:
        initialRfGap = rft.Drift(INITIAL_L)
        if SOLENOID_TYPE == 'HomogeneousChannel':
            initialRfGap.set_static_Bfield(0, 0, SOL_HOMOG_BZ)
        zFinalInVolume += INITIAL_L
        lat.append(initialRfGap)
    tRf = RF_T0
    rfGap = rft.Drift(RF_SEPARATION)
    rfGap.set_aperture(RF_R_APERTURE, RF_R_APERTURE, 'circular')
    if SOLENOID_TYPE == 'HomogeneousChannel':
        rf.set_static_Bfield(0, 0, SOL_HOMOG_BZ)
        rfGap.set_static_Bfield(0, 0, SOL_HOMOG_BZ)
    for structInd in np.arange(RF_N_STRUCTURES):
        print('Setting rf.t0 = {:f} mm/c'.format(tRf))
        rf.set_t0(tRf)
        if structInd == 0:
            rf.set_phid(RF_PHASE_1ST)
            rf.set_P_actual((RF_GRADIENT_1ST/rfttools.RF_CLIC_GRADIENT)**2.)
        else:
            rf.set_phid(RF_PHASE_AFTER_1ST)
            rf.set_P_actual((RF_GRADIENT_AFTER_1ST/rfttools.RF_CLIC_GRADIENT)**2.)
        for rfPeriodInd in np.arange(RF_N_PERIODS_PER_STRUCTURE):
            lat.append(rf)
            # TODO: Get value form RF-Track object rf
        zFinalInVolume += RF_L_STRUCTURE / 2.
        if SOLENOID_TYPE == 'Analytical':
            solZ, solBzOnAxis = bd.generate_solenoid_fieldmap_wilson(
                np.arange(-SOL_L*5., SOL_L*5., 1e-3), 0.,
                SOL_R_IN_COIL, SOL_R_OUT_COIL, SOL_L/2., SOL_J
            )
            solFieldmapStep = solZ[1] - solZ[0]
            solenoid = rft.Static_Magnetic_FieldMap_1d(solBzOnAxis, solFieldmapStep)
            vol.add(solenoid, 0., 0., zFinalInVolume, 'center')
        zFinalInVolume += RF_L_STRUCTURE / 2.
        if splitTracking and structInd == N_RF_STRUCT_1ST_TRACKING - 1:
            zStop1stTracking = zFinalInVolume
        if structInd < RF_N_STRUCTURES-1:
            lat.append(rfGap)
            # TODO: Get value form RF-Track object
            zFinalInVolume += RF_SEPARATION
            tRf += RF_SEPARATION * 1e3   # [mm/c]
    vol.add(lat, 0, 0, amdExitZInVolume, 'entrance')

if FINAL_L > 0:
    finalDrift = rft.Drift(FINAL_L)
    if SOLENOID_TYPE == 'HomogeneousChannel':
        finalDrift.set_static_Bfield(0, 0, SOL_HOMOG_BZ)
    vol.add(finalDrift, 0, 0, zFinalInVolume, 'entrance')
    zFinalInVolume += FINAL_L
if not splitTracking:
    zStop1stTracking = zFinalInVolume

# % Get total field for plotting
# if (1)
#   Z_pl  = linspace(V.get_s0()*1e3, V.get_s1()*1e3, 300); % mm
#   Bz_pl = [];
#   for z = Z_pl
#     [E_pl, B_pl] = V.get_field(0, 0, z, 0);
#     Bz_pl = [ Bz_pl B_pl(3) ];
#   end
#   if (Bz_pl(end)==0)
#     Bz_pl(end) = amd.Bc;
#   endif
#   Z_pl += zv.target_exit + amd.zfte; % in field map coordinate
#   amd.Z_pl = Z_pl;
#   amd.Bz_pl = Bz_pl;
# endif

# % Track with Volume

T0 = rft.TrackingOptions()
T0.dt_mm = 0.2
T0.tt_dt_mm = 1.   # [mm/c], track the emittance every tt_dt_mm (time)
T0.wp_dt_mm = 0.5e3   # [mm/c], save the distr. on disk every wp_dt_mm (time)
# Start tracking at s0
T0.backtrack_at_entrance = False
T0.odeint_algorithm = 'rkf45'   # Options: 'rk2', 'rkf45', 'rk8pd'
T0.odeint_epsabs = 1e-5
T0.open_boundaries = False
T0.verbosity = 1   # 0 (default), 1 or 2

vol.set_s0(0)
vol.set_s1(zStop1stTracking)
T0.t_max_mm = zStop1stTracking * 1e3 + T_ADD_NON_RELATIVISTIC
print('1st particle tracking ends at s1 = {:f} m or at t_max = {:f} mm/c.'.format(
    zStop1stTracking, T0.t_max_mm)
)
B1_6dT = vol.track(B0_6dT, T0)
M1_6dT = B1_6dT.get_phase_space()
B1_6d = vol.get_bunch_at_s1()
M1_6d = B1_6d.get_phase_space("%x %xp %y %yp %t %Pc")
bd.convert_rftrack_to_standard_df(
    rftrackDf=M1_6d, rftrackDfFormat='rftrack_xp_t', pdgId=-11,
    outFwfPath='./Results_AMD/LatestSim/DistrOut_After1stTracking_6d'
)
Bend_6dT = B1_6dT

if splitTracking:
    tCut = np.min(M1_6d[:, 4]) + N_RF_CELLS_LONG_PS_CUT * rfttools.RF_CLIC_L_CELL*1e3  # [mm/c]
    M1_6d_frontBuckets = M1_6d[M1_6d[:, 4] < tCut, :]
    B1_6d.set_phase_space(M1_6d_frontBuckets)
    B1_6dT = rft.Bunch6dT(B1_6d)
    bd.convert_rftrack_to_standard_df(
        rftrackDf=M1_6d_frontBuckets, rftrackDfFormat='rftrack_xp_t', pdgId=-11,
        outFwfPath='./Results_AMD/LatestSim/DistrOut_FrontBuckets_After1stTracking_6d'
    )
    vol.set_s1(zFinalInVolume)
    T0.t_max_mm = zFinalInVolume * 1e3 + T_ADD_NON_RELATIVISTIC
    print('2nd particle tracking ends at s1 = {:f} m or at t_max = {:f} mm/c.'.format(
        zFinalInVolume, T0.t_max_mm)
    )
    B2_6dT = vol.track(B1_6dT, T0)
    B2_6d = vol.get_bunch_at_s1()
    M2_6d = B2_6d.get_phase_space("%x %xp %y %yp %t %Pc")
    bd.convert_rftrack_to_standard_df(
        rftrackDf=M2_6d, rftrackDfFormat='rftrack_xp_t', pdgId=-11,
        outFwfPath='./Results_AMD/LatestSim/DistrOut_After2ndTracking_6d'
    )
    Bend_6dT = B2_6dT
# A_AMD_LOSS = B_AMD_6dT.get_lost_particles();

# TODO: Qbunch
# bd.convert_rftrack_to_standard_df(
#     rftrackDf=M1_6dT, rftrackDfFormat='rftrack_Px_S', pdgId=-11,
#     outFwfPath='./Results_AMD/LatestSim/DistrOut_AMD_6dT'
# )

fig1, ax1 = plt.subplots(5, 1)
rfttools.save_plot_transport(ax1, vol, B0_6dT, Bend_6dT, OUT_REL_PATH)
plt.show(block=False)
input("Press Enter to continue...")
