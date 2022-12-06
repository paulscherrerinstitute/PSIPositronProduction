import RF_Track as rft
import RFTrackTools as rfttools
import BeamDynamics as bd
import SimulationData as sd
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json


# INPUT reproducing YonkeTool_V2, CLIC TW L-band, 0.5 T ###########################################
# BUNCH_FILEPATH = 'RFTrack/YongkeTool_V2/amd_input/E6GeV_SpotSize0.5mm_Target5X0.dat'
# RFTRACK_FORMAT = 'rftrack_xp_t'
# BUNCH_PDGID = -11
# PARTICLE_CHARGE = +1   # [e], +1 = positrons, -1 = electrons
# PARTICLE_MASS = rft.electronmass   # [MeV/c/c]
# BUNCH_Z = 0.
# BUNCH_DOWNSAMPLING = 10
# #
# VOL_R_APERTURE = 1.  # [m]
# #
# TARGET_L = 17.5  # [mm]
# TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = +41e-3   # [m], YonkeTool_V2, CLIC TW L-band, 0.5 T
# #
# AMD_FIELDMAP_2D = 'RFTrack/YongkeTool_V2/field/field_map_HTS_5coils_Apr2022.dat'
# AMD_FIELDMAP_1D = 'RFTrack/YongkeTool_V2/field/field_map_HTS_5coils_Apr2022_1D.dat'
# USE_AMD_FIELDMAP_3D = True
# AMD_R_APERTURE = 20e-3   # [m]
# AMD_L_HALF_MECHANICAL = 96.5e-3   # [m]
# #
# TRACK_AFTER_AMD = True
# #
# INITIAL_L = 0.   # [m]
# #
# # TODO: Use following parameter
# RF_FIELDMAP = None
# RF_FIELDMAP_DIM = '3D'
# RF_FIELDMAP_TYPE = 'SinglePeriod'
# RF_FIELDMAP_GRAD = rfttools.RF_CLIC_GRADIENT
# RF_N_STRUCTURES = 11
# RF_N_CELLS = 30
# RF_N_PERIODS_PER_STRUCTURE = np.round(RF_N_CELLS/rfttools.RF_CLIC_CELLS_PER_PERIOD)
# RF_L_STRUCTURE = RF_N_PERIODS_PER_STRUCTURE * bd.C/rfttools.RF_CLIC_FREQ  # [m]
# RF_L_FLANGE = 0.05461  # [m]
# RF_L_MECH_MARGIN = 0.01  # [m]
# RF_T0 = TARGET_L + 225.3  # [mm/c], optimal for current AMD field map
# RF_PHASES = (171., 171.)  # [deg]
# RF_GRADIENTS = (17.5e6, 21e6)  # [V/m]
# RF_SEPARATION = 0.2   # [m]
# RF_R_APERTURE = 20e-3   # [m]
# #
# AUTOPHASING = False
# #
# SOLENOID_TYPE = 'HomogeneousChannel'
# SOL_HOMOG_BZ = 0.5   # [T]
# # or
# # SOLENOID_TYPE = 'Analytical'
# # SOL_R_IN_COIL = 0.130   # [m]
# # SOL_R_OUT_COIL = 0.250   # [m]
# # SOL_L = RF_L_STRUCTURE - RF_L_FLANGE - rfttools.RF_CLIC_L_CELL - RF_L_MECH_MARGIN
# # SOL_J = 3.54e6   # [A/m2]
# # SOL_HOMOG_BZ = 0.5   # [T]
# #
# FINAL_L = 0.   # [m]
# #
# T_ADD_NON_RELATIVISTIC = 1000.  # [mm/c]
# N_RF_STRUCT_1ST_TRACKING = RF_N_STRUCTURES
# N_RF_CELLS_LONG_PS_CUT = 1 + 2 * rfttools.RF_CLIC_CELLS_PER_PERIOD
# # E.g.: 1 + 2 * ... = keep 3 positron buckets
###################################################################################################


# INPUT reproducing YonkeTool_V3, LargeR TW L-band, 0.5 T #########################################
TRACK_ONLY_REF_PART = False
BUNCH_FILEPATH = 'RFTrack/YongkeTool_V3/Dat/' + \
    'TargetOutputPositrons_E6GeV_SpotSize0.5mm_EmittXY15um_ConvTarget5X0.dat'
RFTRACK_FORMAT = 'rftrack_xp_t'
BUNCH_PDGID = -11
PARTICLE_CHARGE = +1   # [e], +1 = positrons, -1 = electrons
PARTICLE_MASS = rft.electronmass   # [MeV/c/c]
BUNCH_Z = 0.
BUNCH_DOWNSAMPLING = 1
#
VOL_R_APERTURE = 30e-3  # [m]
#
TARGET_L = 17.5  # [mm]
TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = +35e-3   # [m], YonkeTool_V3, LargeR TW L-band, 0.5 T
#
AMD_FIELDMAP_2P5D = 'RFTrack/YongkeTool_V2/field/field_map_HTS_5coils_Apr2022.dat'
AMD_FIELDMAP_1D = 'RFTrack/YongkeTool_V2/field/field_map_HTS_5coils_Apr2022_1D.dat'
USE_AMD_FIELDMAP_2P5D = True
AMD_R_APERTURE = None   # [m]
# TODO: Old value in the following line
AMD_L_HALF_MECHANICAL = 96.5e-3   # [m]
#
TRACK_AFTER_AMD = True
#
ACCEL_WITH_HOMOG_EZ = True
RF_FIELDMAP = 'RFTrack/YongkeTool_V3/field/field_map_LargeR_Lband.dat'
RF_FIELDMAP_DIM = '1D'
RF_FIELDMAP_TYPE = 'Full'
RF_FIELDMAP_GRAD = 20e6  # [V/m]
RF_N_STRUCTURES = 5 + 23  # 5 + 23
RF_L_STRUCTURE = 3.240  # [m]
#   Current RF_L_STRUCTURE including RF_SEPARATION = 3.207 m
RF_FREQ = 2e9  # [Hz]
RF_L_CELL = bd.C / RF_FREQ * 9./20.  # [m]
# RF_L_FLANGE = xxx  # [m]
# RF_L_MECH_MARGIN = xxx  # [m]
# RF_PHASES = (-125.7, -127.8, -132.0, -102.9, -95.0)  # [deg],
#   values for TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = +30e-3 m
RF_PHASES = (-130., -130., -135., -75, -75)  # [deg],
#   values for TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = +35e-3 m
RF_SET_GRADIENTS = (20e6, 20e6, 20e6, 20e6, 20e6)  # [V/m]
RF_SEPARATION = RF_L_STRUCTURE - 3.207140  # [m], struct. separation partially included in fieldmap
RF_R_APERTURE = None   # [m]
#
AUTOPHASING = True
P0_REF = 100.  # [MeV/c]
P0_REF_2 = 205.  # [MeV/c]
T0_REF_2 = 54.92 * bd.C / 1e6  # [mm/c]
#
# SOLENOID_TYPE = 'HomogeneousChannel'
# SOL_HOMOG_BZ = 0.5   # [T]
# or
# SOLENOID_TYPE = 'Analytical'
# SOL_R_IN_COIL = 0.130   # [m]
# SOL_R_OUT_COIL = 0.250   # [m]
# SOL_L = RF_L_STRUCTURE - RF_L_FLANGE - rfttools.RF_CLIC_L_CELL - RF_L_MECH_MARGIN
# SOL_J = 3.54e6   # [A/m2]
# SOL_HOMOG_BZ = 0.5   # [T]
# or
# TODO: 'ZFirstCenter' must usually be changed together with TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD
# SOLENOID_TYPE = 'Simulated'
# SOLENOIDS = {
#     'Type1': {
#         'Type': 'Simulated',
#         'Fieldmap': 'fieldmaps/20220819_Solenoid_type_1_200A.dat',
#         'FieldmapCurrent': 200.,  # [A]
#         'SetCurrent': 200.,  # [A]
#         'MechanicalLength': 200e-3,  # [m]
#         'ZFirstCenter': 565e-3,  # [m]
#         'RepDistances': [321.5e-3, 334e-3],  # [m]
#         'RepNums': [8, 2]
#     },
# }
# SOL_HOMOG_BZ = 0.5  # [T]
# or
SOLENOID_TYPE = 'Simulated'
SOLENOIDS = {
    'Type1a': {
        'Type': 'Simulated',
        'Fieldmap': 'fieldmaps/20220819_Solenoid_type_1_200A.dat',
        'FieldmapCurrent': 200.,  # [A]
        'SetCurrent': 200.,  # [A]
        'MechanicalLength': 200e-3,  # [m]
        'ZFirstCenter': 565e-3,  # [m]
        'RepDistance': 3240e-3,  # [m]
        'RepNum': RF_N_STRUCTURES
    },
    'Type1b': {
        'Type': 'Simulated',
        'Fieldmap': 'fieldmaps/20220819_Solenoid_type_1_200A.dat',
        'FieldmapCurrent': 200.,  # [A]
        'SetCurrent': 200.,  # [A]
        'MechanicalLength': 200e-3,  # [m]
        'ZFirstCenter': 3137e-3,  # [m]
        'RepDistance': 3240e-3,  # [m]
        'RepNum': RF_N_STRUCTURES
    },
    'Type1c': {
        'Type': 'Simulated',
        'Fieldmap': 'fieldmaps/20220819_Solenoid_type_1_200A.dat',
        'FieldmapCurrent': 200.,  # [A]
        'SetCurrent': 200.,  # [A]
        'MechanicalLength': 200e-3,  # [m]
        'ZFirstCenter': 3471e-3,  # [m]
        'RepDistance': 3240e-3,  # [m]
        'RepNum': RF_N_STRUCTURES
    },
    'Type2': {
        'Type': 'Simulated',
        'Fieldmap': 'fieldmaps/20220819_Solenoid_type_2_200A.dat',
        'FieldmapCurrent': 200.,  # [A]
        'SetCurrent': 200.,  # [A]
        'MechanicalLength': 2214e-3,  # [m]
        'ZFirstCenter': 1851e-3,  # [m]
        'RepDistance': 3240e-3,  # [m]
        'RepNum': RF_N_STRUCTURES
    },
    'Tuning': {
        'Type': 'Simulated',
        'Fieldmap': 'fieldmaps/20221011_Solenoid_Tuning_V1.dat',
        'FieldmapCurrent': 200.,  # [A]
        'SetCurrent': 200.,  # [A]
        'MechanicalLength': 72e-3,  # [m]
        'ZFirstCenter': 300e-3,  # [m]
        'RepDistance': 0,  # [m]
        'RepNum': 1
    }
}
SOL_HOMOG_BZ = 0.5  # [T]
#
INITIAL_L = RF_SEPARATION / 2.  # [m]
FINAL_L = 1.  # [m]
#
T_ADD_NON_RELATIVISTIC = 1000.  # [mm/c]
N_RF_STRUCT_1ST_TRACKING = 5
N_RF_CELLS_LONG_PS_CUT = 2
PZ_MIN_LONG_PS_CUT = 40.  # [MeV/c]
###################################################################################################


splitTracking = RF_N_STRUCTURES > N_RF_STRUCT_1ST_TRACKING

OUT_REL_PATH = './RFTrackOutput/LatestSimCaptureLinac/'

# TODO: Refactorize, same code in Run_Linac1_Section1_Simple.py
distrMatNp = np.loadtxt(BUNCH_FILEPATH, skiprows=1)
beamIn, _ = bd.convert_rftrack_to_standard_df(
    rftrackDf=distrMatNp[::BUNCH_DOWNSAMPLING, :], rftrackDfFormat=RFTRACK_FORMAT,
    s=BUNCH_Z, pdgId=BUNCH_PDGID, Qbunch=bd.PART_CONSTS['Q'][BUNCH_PDGID]*distrMatNp.shape[0]
)
M0 = bd.convert_standard_df_to_rftrack(
    standardDf=beamIn, rftrackDfFormat=RFTRACK_FORMAT
)[0].to_numpy()
# TODO: Verify following change
BUNCH_POPULATION = M0.shape[0]  # YongkeTool_V2
# BUNCH_POPULATION = 0  # YongkeTool_V3, but only in track_AMD_HTS.m
B0_6d = rft.Bunch6d(PARTICLE_MASS, BUNCH_POPULATION, PARTICLE_CHARGE, M0)
B0_6dT = rft.Bunch6dT(B0_6d)

TARGET_EXIT_Z_IN_VOLUME = 0.   # [m]
vol = rft.Volume()
vol.set_aperture(VOL_R_APERTURE, VOL_R_APERTURE, 'circular')
beamlineSetup = pd.DataFrame(
    columns=['ElementType', 'zWrtTargetExit', 'MechanicalLength', 'Fieldmap']
)

# TODO: Clean distinction between AMD_FIELDMAP_1D AND _2D
if USE_AMD_FIELDMAP_2P5D:
    amdFieldmap = bd.load_octave_matrices(AMD_FIELDMAP_2P5D)
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
amdFieldLength = amdZ[indZConstBzStart] - amdZ[indZTargetExit]  # [mm]
if SOLENOID_TYPE == 'HomogeneousChannel':
    indsEff = np.arange(indZTargetExit, indZConstBzStart+1)
elif SOLENOID_TYPE in ['Analytical', 'Simulated']:
    indsEff = np.arange(indZTargetExit, amdBzOnAxis.shape[0])
    amdFieldLengthAnalytical = \
        amdZ[amdBzOnAxis.shape[0]-1] - amdZ[indZTargetExit]
amdZEff = amdZ[indsEff]
amdExitZInVolume = TARGET_EXIT_Z_IN_VOLUME + amdFieldLength*1e-3   # [m]
print('amdFieldLength = {:f} mm.'.format(amdFieldLength))
if USE_AMD_FIELDMAP_2P5D:
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
elif SOLENOID_TYPE in ['Analytical', 'Simulated']:
    amd.set_length(amdFieldLengthAnalytical*1e-3)
    # TODO: Is the following line necessary in volume?
    amd.set_nsteps(int(amdFieldLengthAnalytical*2.))
vol.add(amd, 0, 0, TARGET_EXIT_Z_IN_VOLUME, 'entrance')
beamlineSetup.loc[len(beamlineSetup.index)] = [
    'AMD',
    TARGET_EXIT_Z_IN_VOLUME - TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD,
    AMD_L_HALF_MECHANICAL * 2.,
    os.path.basename(AMD_FIELDMAP_2P5D)
]
beamlineSetup.loc[len(beamlineSetup.index)] = [
    'TargetExit', TARGET_EXIT_Z_IN_VOLUME, 0., ''
]
zFinalInVolume = amdExitZInVolume

AMD_MECH_EXIT_Z_IN_VOLUME = \
    TARGET_EXIT_Z_IN_VOLUME - TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD + AMD_L_HALF_MECHANICAL
if AMD_R_APERTURE is not None:
    amdAperture = rft.Drift(AMD_MECH_EXIT_Z_IN_VOLUME)
    amdAperture.set_aperture(AMD_R_APERTURE, AMD_R_APERTURE, 'circular')
    vol.add(amdAperture, 0, 0, 0, 'entrance')

if TRACK_AFTER_AMD:
    lat = rft.Lattice()
    if INITIAL_L > 0:
        initialRfGap = rft.Drift(INITIAL_L)
        if SOLENOID_TYPE == 'HomogeneousChannel':
            initialRfGap.set_static_Bfield(0, 0, SOL_HOMOG_BZ)
        zFinalInVolume += INITIAL_L
        lat.append(initialRfGap)
    if not AUTOPHASING:
        tRf = RF_T0  # [mm/c]
    else:
        t0RefPart1 = TARGET_L + amdFieldLength  # [mm/c]
        refPart1Vol = np.array([
            0., 0., 0., 0., amdFieldLength, P0_REF, PARTICLE_MASS, PARTICLE_CHARGE, +1., t0RefPart1
        ])
        refPart1Lat = np.array([
            0., 0., 0., 0., t0RefPart1, P0_REF, PARTICLE_MASS, PARTICLE_CHARGE, +1.
        ])
        tRf = None
    rfGap = rft.Drift()
    if RF_R_APERTURE is not None:
        rfGap.set_aperture(RF_R_APERTURE, RF_R_APERTURE, 'circular')
    if SOLENOID_TYPE == 'HomogeneousChannel':
        rfGap.set_static_Bfield(0, 0, SOL_HOMOG_BZ)
        rfHomogBz = SOL_HOMOG_BZ
    else:
        rfHomogBz = None
    for structInd in np.arange(RF_N_STRUCTURES):
        try:
            powerScalingFactor = (RF_SET_GRADIENTS[structInd] / RF_FIELDMAP_GRAD) ** 2.
        except IndexError:
            powerScalingFactor = (RF_SET_GRADIENTS[-1] / RF_FIELDMAP_GRAD) ** 2.
        try:
            rfPhase = RF_PHASES[structInd]
            # rfPhase = 0
        except IndexError:
            rfPhase = RF_PHASES[-1]
            # rfPhase = 0
        try:
            print('Setting rf.t0 = {:f} mm/c'.format(tRf))
        except TypeError:
            print('Using autophasing...')
        if RF_FIELDMAP_TYPE == 'SinglePeriod':
            # TODO: Load fieldmaps only once.
            rf = rfttools.rf_struct_from_single_period(
                RF_FIELDMAP, RF_FIELDMAP_DIM, RF_N_PERIODS_PER_STRUCTURE, powerScalingFactor, tRf, rfPhase,
                aperture=RF_R_APERTURE, additionalHomogBz=rfHomogBz
            )
        elif RF_FIELDMAP_TYPE == 'Full':
            rf = rfttools.rf_struct_from_full_fieldmap(
                RF_FIELDMAP, RF_FIELDMAP_DIM, powerScalingFactor, tRf, rfPhase,
                aperture=RF_R_APERTURE, additionalHomogBz=rfHomogBz
            )
        if ACCEL_WITH_HOMOG_EZ and structInd > N_RF_STRUCT_1ST_TRACKING - 1:
            rf = rft.Drift(rf.get_length())
            rf.set_static_Efield(0, 0, RF_SET_GRADIENTS[-1])
        lat.append(rf)
        zFinalInVolume += rf.get_length() / 2.
        beamlineSetup.loc[len(beamlineSetup.index)] = [
            'RF', zFinalInVolume, RF_L_STRUCTURE, os.path.basename(RF_FIELDMAP)
        ]
        if SOLENOID_TYPE == 'Analytical':
            solenoid = rfttools.solenoid_from_analytical_formula(
                SOL_L, SOL_R_IN_COIL, SOL_R_OUT_COIL, SOL_J
            )
            vol.add(solenoid, 0., 0., zFinalInVolume, 'center')
        zFinalInVolume += rf.get_length() / 2.
        if splitTracking and structInd == N_RF_STRUCT_1ST_TRACKING - 1:
            zStop1stTracking = zFinalInVolume
        if structInd < RF_N_STRUCTURES-1:
            rfGap.set_length(RF_SEPARATION)
        else:
            rfGap.set_length(RF_SEPARATION / 2.)
        lat.append(rfGap)
        zFinalInVolume += rfGap.get_length()
        if not AUTOPHASING:
            tRf += rfGap.get_length() * 1e3   # [mm/c]
    if AUTOPHASING:
        print('Autophasing in lattice...')
        pzFinalAutophasing = lat.autophase(rft.Bunch6d(refPart1Lat))
        print('pzFinalAutophasing 1st tracking = {:e} MeV/c'.format(pzFinalAutophasing))
    vol.add(lat, 0, 0, amdExitZInVolume, 'entrance')
    # TODO: Simplify script wrt solenoid (analytical or simulated is almost the same)
    for solType, sol in SOLENOIDS.items():
        if sol['Type'] == 'Simulated':
            solenoidCenters = [sol['ZFirstCenter'], ]
            for repInd in np.arange(sol['RepNum'] - 1):
                solenoidCenters.append(solenoidCenters[-1] + sol['RepDistance'])
            solenoid = rfttools.solenoid_from_fieldmap(
                sol['Fieldmap'], sol['FieldmapCurrent'], sol['SetCurrent']
            )
            for zCenter in solenoidCenters:
                vol.add(solenoid, 0., 0., zCenter, 'center')
                beamlineSetup.loc[len(beamlineSetup.index)] = [
                    solType, zCenter, sol['MechanicalLength'], os.path.basename(sol['Fieldmap'])
                ]
    # if AUTOPHASING:
    #     print('Autophasing in volume...')
    #     pzFinalAutophasing = vol.autophase(rft.Bunch6dT(refPart1Vol))
    #     print('pzFinalAutophasing 1st tracking = {:e} MeV/c'.format(pzFinalAutophasing))

if FINAL_L > 0:
    finalDrift = rft.Drift(FINAL_L)
    if SOLENOID_TYPE == 'HomogeneousChannel':
        finalDrift.set_static_Bfield(0, 0, SOL_HOMOG_BZ)
    vol.add(finalDrift, 0, 0, zFinalInVolume, 'entrance')
    zFinalInVolume += finalDrift.get_length()
if not splitTracking:
    zStop1stTracking = zFinalInVolume

beamlineSetup['zWrtAmdPeakField'] = \
    beamlineSetup['zWrtTargetExit'] + TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD
beamlineSetup[[
    'ElementType', 'zWrtTargetExit', 'zWrtAmdPeakField', 'MechanicalLength', 'Fieldmap'
]].to_csv(os.path.join(OUT_REL_PATH, 'BeamlineSetup.dat'), index=None)

trackingOpts = rft.TrackingOptions()
trackingOpts.dt_mm = 0.2
trackingOpts.tt_dt_mm = 1.   # [mm/c], track the emittance every tt_dt_mm (time)
trackingOpts.wp_dt_mm = 0.5e3   # [mm/c], save the distr. on disk every wp_dt_mm (time)
trackingOpts.backtrack_at_entrance = False
trackingOpts.odeint_algorithm = 'rkf45'   # Options: 'rk2', 'rkf45', 'rk8pd'
trackingOpts.odeint_epsabs = 1e-5
trackingOpts.verbosity = 1   # 0 (default), 1 or 2

vol.set_s0(BUNCH_Z)
vol.set_s1(zStop1stTracking)
trackingOpts.t_max_mm = zStop1stTracking * 1e3 + T_ADD_NON_RELATIVISTIC
print('1st particle tracking ends at s1 = {:f} m or at t_max = {:f} mm/c.'.format(
    zStop1stTracking, trackingOpts.t_max_mm)
)
if TRACK_ONLY_REF_PART:
    B1_6dT = vol.track(rft.Bunch6dT(refPart1Vol), trackingOpts)
else:
    B1_6dT = vol.track(B0_6dT, trackingOpts)
M1_6dT = B1_6dT.get_phase_space()
B1_6d = vol.get_bunch_at_s1()
M1_6d = B1_6d.get_phase_space("%x %xp %y %yp %t %Pc")
bd.convert_rftrack_to_standard_df(
    rftrackDf=M1_6d, rftrackDfFormat='rftrack_xp_t', s=B1_6d.get_S()*1e3, pdgId=BUNCH_PDGID,
    outFwfPath=os.path.join(OUT_REL_PATH, 'DistrOut_After1stTracking_6d')
)

fig1, ax1 = plt.subplots(7, 1)
rfttools.save_plot_transport(ax1, vol, B0_6dT, B1_6dT, OUT_REL_PATH, outSuffix='1')
plt.show(block=False)

if splitTracking:
    refPart2 = np.array([
        0., 0., 0., 0., zStop1stTracking*1e3, P0_REF_2,
        PARTICLE_MASS, PARTICLE_CHARGE, +1., T0_REF_2
    ])
    # if AUTOPHASING:
    #     pzFinalAutophasing = vol.autophase(rft.Bunch6dT(refPart2))
    #     print('pzFinalAutophasing 2nd tracking = {:e} MeV/c'.format(pzFinalAutophasing))
    tCut = np.min(M1_6d[:, 4]) + N_RF_CELLS_LONG_PS_CUT * RF_L_CELL*1e3  # [mm/c]
    M1_6d_frontBuckets = M1_6d[(M1_6d[:, 4] < tCut) & (M1_6d[:, 5] > PZ_MIN_LONG_PS_CUT), :]
    B1_6d.set_phase_space(M1_6d_frontBuckets)
    # B1_6dT = rft.Bunch6dT(refPart2)
    # or
    B1_6dT = rft.Bunch6dT(B1_6d)
    bd.convert_rftrack_to_standard_df(
        rftrackDf=M1_6d_frontBuckets, rftrackDfFormat='rftrack_xp_t',
        s=B1_6d.get_S()*1e3, pdgId=BUNCH_PDGID,
        outFwfPath=os.path.join(OUT_REL_PATH, 'DistrOut_FrontBuckets_After1stTracking_6d')
    )
    vol.set_s1(zFinalInVolume)
    trackingOpts.t_max_mm = zFinalInVolume * 1e3 + T_ADD_NON_RELATIVISTIC
    print('2nd particle tracking ends at s1 = {:f} m or at t_max = {:f} mm/c.'.format(
        zFinalInVolume, trackingOpts.t_max_mm)
    )
    B2_6dT = vol.track(B1_6dT, trackingOpts)
    B2_6d = vol.get_bunch_at_s1()
    M2_6d = B2_6d.get_phase_space("%x %xp %y %yp %t %Pc")
    bd.convert_rftrack_to_standard_df(
        rftrackDf=M2_6d, rftrackDfFormat='rftrack_xp_t', s=B2_6d.get_S()*1e3, pdgId=BUNCH_PDGID,
        outFwfPath=os.path.join(OUT_REL_PATH, 'DistrOut_After2ndTracking_6d')
    )

    rfttools.save_plot_transport(ax1, vol, B1_6dT, B2_6dT, OUT_REL_PATH, outSuffix='2')
    plt.show(block=False)

input("Press Enter to continue...")
