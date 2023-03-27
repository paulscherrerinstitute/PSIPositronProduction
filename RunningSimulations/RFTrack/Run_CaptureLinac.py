import RF_Track as rft
import RFTrackTools as rfttools
import OctavePythonInterface as opi
import BeamDynamics as bd
import SimulationData as sd
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json


# INPUT reproducing YonkeTool_V2, CLIC TW L-band, 0.5 T ###########################################
# TRACK_ONLY_REF_PART = False
# BUNCH_FILEPATH = 'RunningSimulations/RFTrack/YongkeTool_V2/amd_input/' + \
#     'E6GeV_SpotSize0.5mm_Target5X0.dat'
# RFTRACK_FORMAT = 'rftrack_xp_t'
# BUNCH_PDGID = -11
# PARTICLE_CHARGE = +1  # [e], +1 = positrons, -1 = electrons
# PARTICLE_MASS = rft.electronmass   # [MeV/c/c]
# BUNCH_Z = 0.
# BUNCH_DOWNSAMPLING = 20
# #
# VOL_R_APERTURE = 1.  # [m]
# #
# TARGET_L = 17.5  # [mm]
# TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = +41e-3  # [m], YonkeTool_V2, CLIC TW L-band, 0.5 T
# #
# AMD_FIELDMAP_2P5D = 'RunningSimulations/RFTrack/YongkeTool_V2/field/' + \
#     'field_map_HTS_5coils_Apr2022.dat'
# AMD_FIELDMAP_1D = 'RunningSimulations/RFTrack/YongkeTool_V2/field/' + \
#     'field_map_HTS_5coils_Apr2022_1D.dat'
# USE_AMD_FIELDMAP_2P5D = True
# AMD_R_APERTURE = 20e-3   # [m]
# AMD_L_HALF_MECHANICAL = 96.5e-3   # [m]
# #
# TRACK_AFTER_AMD = True
# #
# ACCEL_WITH_HOMOG_EZ = False
# RF_FIELDMAP = 'RunningSimulations/RFTrack/YongkeTool_V2/field/field_map_CLIC_Lband.dat'
# RF_FIELDMAP_DIM = '3D_CylindricalSym'
# RF_FIELDMAP_TYPE = 'SinglePeriod'
# RF_FIELDMAP_GRAD = 11.23e6  # [V/m]
# RF_N_STRUCTURES = 11  # 11
# RF_N_CELLS = 30
# RF_CLIC_CELLS_PER_PERIOD = 3.
# RF_CLIC_FREQ = 1.9986163867e+09  # [Hz]
# RF_N_PERIODS_PER_STRUCTURE = np.round(RF_N_CELLS/RF_CLIC_CELLS_PER_PERIOD)
# RF_L_STRUCTURE = RF_N_PERIODS_PER_STRUCTURE * bd.C / RF_CLIC_FREQ  # [m]
# RF_L_FLANGE = 0.05461  # [m]
# RF_L_MECH_MARGIN = 0.01  # [m]
# RF_T0 = TARGET_L + 225.3  # [mm/c], optimal for current AMD field map
# RF_PHASES = (171., 171.)  # [deg]
# RF_SET_GRADIENTS = (17.5e6, 21e6)  # [V/m]
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
# INITIAL_L = 0.  # [m]
# FINAL_L = 0.  # [m]
# #
# T_ADD_NON_RELATIVISTIC = 1000.  # [mm/c]
# N_RF_STRUCT_1ST_TRACKING = RF_N_STRUCTURES
# N_RF_CELLS_LONG_PS_CUT = 1 + 2 * RF_CLIC_CELLS_PER_PERIOD
# # E.g.: 1 + 2 * ... = keep 3 positron buckets
###################################################################################################


# INPUT reproducing YonkeTool_V3, LargeR TW L-band, 0.5 T #########################################
TRACK_ONLY_REF_PART = False
BUNCH_FILEPATH = 'RunningSimulations/RFTrack/YongkeTool_V3/Dat/' + \
    'TargetOutputPositrons_E6GeV_SpotSize0.5mm_EmittXY15um_ConvTarget5X0.dat'
RFTRACK_FORMAT = 'rftrack_xp_t'
BUNCH_PDGID = -11
PARTICLE_CHARGE = +1   # [e], +1 = positrons, -1 = electrons
PARTICLE_MASS = rft.electronmass   # [MeV/c/c]
BUNCH_Z = 0.
BUNCH_DOWNSAMPLING = 20
#
VOL_R_APERTURE = None  # [m]
#
TARGET_L = 17.5  # [mm]
TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = +35e-3   # [m], YonkeTool_V3, LargeR TW L-band, 0.5 T
#
AMD_FIELDMAP_2P5D = 'RunningSimulations/RFTrack/YongkeTool_V2/field/' + \
    'field_map_HTS_5coils_Apr2022.dat'
AMD_FIELDMAP_1D = 'RunningSimulations/RFTrack/YongkeTool_V2/field/' + \
    'field_map_HTS_5coils_Apr2022_1D.dat'
USE_AMD_FIELDMAP_2P5D = True
AMD_R_APERTURE = 30e-3   # [m]
# TODO: Old value in the following line
AMD_L_HALF_MECHANICAL = 96.5e-3   # [m]
#
TRACK_AFTER_AMD = True
#
ACCEL_WITH_HOMOG_EZ = False
#
# RF_FIELDMAP = 'RunningSimulations/RFTrack/YongkeTool_V3/field/field_map_LargeR_Lband.dat'
# RF_FIELDMAP_DIM = '1D'
# RF_FIELDMAP_TYPE = 'Full'
# RF_PHASE_CORR = 0  # [deg]
# or
RF_FIELDMAP = 'Data/Fieldmaps/pLinacF3_full44cells_YZplane_dy2mm_dz0p1L.dat'
RF_FIELDMAP_DIM = '2D'
RF_FIELDMAP_TYPE = 'Full'
RF_PHASE_CORR = -0.4  # [deg]
#
RF_FIELDMAP_GRAD = 20e6  # [V/m]
RF_N_STRUCTURES = 28  # 5 + 23
RF_L_STRUCTURE = 3.240  # [m]
#   Current RF_L_STRUCTURE including RF_SEPARATION = 3.207 m
RF_FREQ = 2e9  # [Hz]
RF_L_CELL = bd.C / RF_FREQ * 9./20.  # [m]
# RF_L_FLANGE = xxx  # [m]
# RF_L_MECH_MARGIN = xxx  # [m]
# RF_PHASES = (-125.7, -127.8, -132.0, -102.9, -95.0)  # [deg],
#   values for TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = +30e-3 m
RF_PHASES = np.array([-130., -130., -135., -75, -75]) + RF_PHASE_CORR  # [deg]
#   values for TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = +35e-3 m
RF_SET_GRADIENTS = (20e6, 20e6, 20e6, 20e6, 20e6)  # [V/m]
RF_R_APERTURE = 30e-3  # [m]
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
#         'Fieldmap': 'Data/Fieldmaps/20220819_Solenoid_type_1_200A.dat',
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
        'Fieldmap': 'Data/Fieldmaps/20220819_Solenoid_type_1_200A.dat',
        'FieldmapCurrent': 200.,  # [A]
        'SetCurrent': 200.,  # [A]
        'MechanicalLength': 200e-3,  # [m]
        'ZFirstCenter': 565e-3,  # [m]
        'RepDistance': 3240e-3,  # [m]
        'RepNum': RF_N_STRUCTURES
    },
    'Type1b': {
        'Type': 'Simulated',
        'Fieldmap': 'Data/Fieldmaps/20220819_Solenoid_type_1_200A.dat',
        'FieldmapCurrent': 200.,  # [A]
        'SetCurrent': 200.,  # [A]
        'MechanicalLength': 200e-3,  # [m]
        'ZFirstCenter': 3137e-3,  # [m]
        'RepDistance': 3240e-3,  # [m]
        'RepNum': RF_N_STRUCTURES
    },
    'Type1c': {
        'Type': 'Simulated',
        'Fieldmap': 'Data/Fieldmaps/20220819_Solenoid_type_1_200A.dat',
        'FieldmapCurrent': 200.,  # [A]
        'SetCurrent': 200.,  # [A]
        'MechanicalLength': 200e-3,  # [m]
        'ZFirstCenter': 3471e-3,  # [m]
        'RepDistance': 3240e-3,  # [m]
        'RepNum': RF_N_STRUCTURES
    },
    'Type2': {
        'Type': 'Simulated',
        'Fieldmap': 'Data/Fieldmaps/20220819_Solenoid_type_2_200A.dat',
        'FieldmapCurrent': 200.,  # [A]
        'SetCurrent': 200.,  # [A]
        'MechanicalLength': 2214e-3,  # [m]
        'ZFirstCenter': 1851e-3,  # [m]
        'RepDistance': 3240e-3,  # [m]
        'RepNum': RF_N_STRUCTURES
    },
    'Tuning': {
        'Type': 'Simulated',
        'Fieldmap': 'Data/Fieldmaps/20221011_Solenoid_Tuning_V1.dat',
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
CHICANE_INSERT = True
CHICANE_FIELDMAP = 'Data/Fieldmaps/field_map_chicane_all.dat'
CHICANE_FIELD_PEAK = 0.1  # [T]
CHICANE_AFTER_RF_STRUCT_NO = 5
CHICANE_TOT_LENGTH = 3.0  # [m]
CHICANE_BEAM_PIPE_HALF_APERTURE_X = 0.075  # [m]
CHICANE_BEAM_PIPE_HALF_APERTURE_Y = 0.020  # [m]
CHICANE_COLLIM_X_INSERT = False
CHICANE_COLLIM_X_LENGTH = 0.12  # [m]
CHICANE_COLLIM_X_TOT_APERTURE = 0.050  # [m]
CHICANE_COLLIM_X_OFFSET = -0.025  # [m]
CHICANE_COLLIM_X_Z_FROM_CENTER = 0.1325  # [m]
#
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
BUNCH_POPULATION = M0.shape[0]
B0_6d = rft.Bunch6d(PARTICLE_MASS, BUNCH_POPULATION, PARTICLE_CHARGE, M0)
B0_6dT = rft.Bunch6dT(B0_6d)

TARGET_EXIT_Z_IN_VOLUME = 0.   # [m]
vol = rft.Volume()
if VOL_R_APERTURE is not None:
    vol.set_aperture(VOL_R_APERTURE, VOL_R_APERTURE, 'circular')
beamlineSetup = pd.DataFrame(
    columns=['ElementType', 'zWrtTargetExit', 'MechanicalLength', 'Fieldmap']
)

# TODO: Clean distinction between AMD_FIELDMAP_1D AND _2D
if USE_AMD_FIELDMAP_2P5D:
    amdFieldmap = opi.load_octave_matrices(AMD_FIELDMAP_2P5D)
    amdDz = amdFieldmap['Z'][0, 1] - amdFieldmap['Z'][0, 0]  # [mm]
    amdDr = amdFieldmap['R'][1, 0]-amdFieldmap['R'][0, 0]  # [mm]
    amdZ = amdFieldmap['Z'][0, :]  # [mm]
    amdBzOnAxis = amdFieldmap['Bz'][0, :]
else:
    amdFieldmap = opi.load_octave_matrices(AMD_FIELDMAP_1D)
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
    # TODO: Very probably redundant
    amd.set_length(amdFieldLength*1e-3)
elif SOLENOID_TYPE in ['Analytical', 'Simulated']:
    amd.set_length(amdFieldLengthAnalytical*1e-3)
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

rfField = opi.load_octave_matrices(RF_FIELDMAP)
rfSeparation = RF_L_STRUCTURE - np.diff(rfField['Z'].flatten()[[0, -1]])[0]
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
            rf = rfttools.rf_struct_from_single_period(
                rfField, RF_FIELDMAP_DIM, RF_N_PERIODS_PER_STRUCTURE, powerScalingFactor,
                tRf, rfPhase, aperture=RF_R_APERTURE, additionalHomogBz=rfHomogBz)
        elif RF_FIELDMAP_TYPE == 'Full':
            rf = rfttools.rf_from_field_map(
                rfField, RF_FIELDMAP_DIM, powerScalingFactor,
                tRf, rfPhase, aperture=RF_R_APERTURE, additionalHomogBz=rfHomogBz)
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
        if structInd == RF_N_STRUCTURES-1 or (
                CHICANE_INSERT and
                structInd in [CHICANE_AFTER_RF_STRUCT_NO-1, CHICANE_AFTER_RF_STRUCT_NO]):
            rfGap.set_length(rfSeparation / 2.)
        else:
            rfGap.set_length(rfSeparation)
        lat.append(rfGap)
        zFinalInVolume += rfGap.get_length()
        if not AUTOPHASING:
            tRf += rfGap.get_length() * 1e3  # [mm/c]
        if CHICANE_INSERT and structInd == CHICANE_AFTER_RF_STRUCT_NO-1:
            # TODO: Refactor insertion of gap
            chicaneGap = rft.Drift(CHICANE_TOT_LENGTH)
            lat.append(chicaneGap)
            chicaneCenter = zFinalInVolume + chicaneGap.get_length()/2.
            zFinalInVolume += chicaneGap.get_length()
            if not AUTOPHASING:
                tRf += chicaneGap.get_length() * 1e3  # [mm/c]
    if AUTOPHASING:
        print('Autophasing in lattice...')
        pzFinalAutophasing = lat.autophase(rft.Bunch6d(refPart1Lat))
        print('pzFinalAutophasing 1st tracking = {:e} MeV/c'.format(pzFinalAutophasing))
    vol.add(lat, 0, 0, amdExitZInVolume, 'entrance')
    # TODO: Simplify script wrt solenoid (analytical or simulated is almost the same)
    if SOLENOID_TYPE != 'HomogeneousChannel':
        for solType, sol in SOLENOIDS.items():
            if sol['Type'] == 'Simulated':
                solenoidCenters = [sol['ZFirstCenter'], ]
                for repInd in np.arange(sol['RepNum'] - 1):
                    solenoidCenters.append(solenoidCenters[-1] + sol['RepDistance'])
                solenoidCenters = np.array(solenoidCenters)
                if CHICANE_INSERT:
                    solenoidCenters[CHICANE_AFTER_RF_STRUCT_NO:] += CHICANE_TOT_LENGTH
                solenoid = rfttools.solenoid_from_fieldmap(
                    sol['Fieldmap'], sol['FieldmapCurrent'], sol['SetCurrent']
                )
                for zCenter in solenoidCenters:
                    vol.add(solenoid, 0., 0., zCenter, 'center')
                    beamlineSetup.loc[len(beamlineSetup.index)] = [
                        solType, zCenter, sol['MechanicalLength'], os.path.basename(sol['Fieldmap'])
                    ]
    if CHICANE_INSERT:
        chicaneField = opi.load_octave_matrices(CHICANE_FIELDMAP)
        strengthFactor = CHICANE_FIELD_PEAK / chicaneField['field_peak']
        dx = chicaneField['X'][1, 0, 0] - chicaneField['X'][0, 0, 0]
        dy = chicaneField['Y'][0, 1, 0] - chicaneField['Y'][0, 0, 0]
        dz = chicaneField['Z'][0, 0, 1] - chicaneField['Z'][0, 0, 0]
        magnetL = chicaneField['Z'][0, 0, -1] - chicaneField['Z'][0, 0, 0]
        # Variant 1: rft.Static_Magnetic_FieldMap() ensures div(B) = 0
        chicaneDipoles = rft.Static_Magnetic_FieldMap(
            chicaneField['Bx']*strengthFactor, chicaneField['By']*strengthFactor,
            chicaneField['Bz']*strengthFactor, chicaneField['X'][0, 0, 0],
            chicaneField['Y'][0, 0, 0], dx, dy, dz, magnetL)
        # Variant 2: rft.RF_FieldMap()
        # vanishingE = np.zeros(chicaneField['Bx'].shape)
        # chicaneDipoles = rft.RF_FieldMap(
        #     vanishingE, vanishingE, vanishingE, chicaneField['Bx']*strengthFactor*1j,
        #     chicaneField['By']*strengthFactor*1j, chicaneField['Bz']*strengthFactor*1j,
        #     chicaneField['X'][0, 0, 0], chicaneField['Y'][0, 0, 0], dx, dy, dz, magnetL, 0., 1)
        # chicaneDipoles.set_t0(0)
        # chicaneDipoles.set_phid(-90.)
        vol.add(chicaneDipoles, 0, 0, chicaneCenter, 'center')
        chicaneBeamPipe = rft.Drift(CHICANE_TOT_LENGTH)
        chicaneBeamPipe.set_aperture(
            CHICANE_BEAM_PIPE_HALF_APERTURE_X, CHICANE_BEAM_PIPE_HALF_APERTURE_Y, 'rectangular')
        vol.add(chicaneBeamPipe, 0, 0, chicaneCenter, 'center')
        if CHICANE_COLLIM_X_INSERT:
            chicaneCollimX = rft.Drift(CHICANE_COLLIM_X_LENGTH)
            chicaneCollimX.set_aperture(CHICANE_COLLIM_X_TOT_APERTURE/2., np.Inf, 'rectangular')
            vol.add(
                chicaneCollimX, CHICANE_COLLIM_X_OFFSET, 0,
                chicaneCenter+CHICANE_COLLIM_X_Z_FROM_CENTER, 'entrance')
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

fig1, ax1 = plt.subplots(9, 1)
rfttools.save_plot_transport(ax1, vol, B0_6dT, B1_6dT, OUT_REL_PATH, outSuffix='1')
plt.show(block=False)

yMesh = np.arange(0, 0.03, 1e-3)
zMesh = np.arange(0.25, 1., 1e-3)
emFields = rfttools.save_em_fields(vol, [0], yMesh, zMesh, returnMultidimNpArray=True)
r = emFields[0, :, :, 1]
z = emFields[0, :, :, 2]
Er = emFields[0, :, :, 4]
Etheta = emFields[0, :, :, 3]
Ez = emFields[0, :, :, 5]
Br = emFields[0, :, :, 7]
Btheta = emFields[0, :, :, 6]
Bz = emFields[0, :, :, 8]
BzSolenoids = 0.5  # [T], attention, this is only the homogeneous part!
Emax = 35e6  # [ V/m]

fig2, ax2 = plt.subplots(3, 2)
p00 = ax2[0, 0].pcolormesh(z, r, Er, shading='nearest', vmin=-Emax, vmax=Emax)
c00 = fig2.colorbar(p00, ax=ax2[0, 0])
ax2[0, 0].set_xlabel('z [m]')
ax2[0, 0].set_ylabel('r [m]')
c00.set_label('Er [V/m]', rotation=270)
p10 = ax2[1, 0].pcolormesh(z, r, Etheta, shading='nearest', vmin=-Emax, vmax=Emax)
c10 = fig2.colorbar(p10, ax=ax2[1, 0])
ax2[1, 0].set_xlabel('z [m]')
ax2[1, 0].set_ylabel('r [m]')
c10.set_label('Etheta [V/m]', rotation=270)
p20 = ax2[2, 0].pcolormesh(z, r, Ez, shading='nearest', vmin=-Emax, vmax=Emax)
c20 = fig2.colorbar(p20, ax=ax2[2, 0])
ax2[2, 0].set_xlabel('z [m]')
ax2[2, 0].set_ylabel('r [m]')
c20.set_label('Ez [V/m]', rotation=270)
p01 = ax2[0, 1].pcolormesh(z, r, Br, shading='nearest', vmin=-0.07, vmax=0.07)
c01 = fig2.colorbar(p01, ax=ax2[0, 1])
ax2[0, 1].set_xlabel('z [m]')
ax2[0, 1].set_ylabel('r [m]')
c01.set_label('Br [T]', rotation=270)
p11 = ax2[1, 1].pcolormesh(z, r, Btheta, shading='nearest', vmin=-0.07, vmax=0.07)
c11 = fig2.colorbar(p11, ax=ax2[1, 1])
ax2[1, 1].set_xlabel('z [m]')
ax2[1, 1].set_ylabel('r [m]')
c11.set_label('Btheta [T]', rotation=270)
p21 = ax2[2, 1].pcolormesh(z, r, Bz-BzSolenoids, shading='nearest', vmin=-0.07, vmax=0.07)
c21 = fig2.colorbar(p21, ax=ax2[2, 1])
ax2[2, 1].set_xlabel('z [m]')
ax2[2, 1].set_ylabel('r [m]')
c21.set_label('Bz [T]', rotation=270)
plt.show(block=False)

rInds = [0, 15, 29]
fig3, ax3 = plt.subplots(2, 1)
for rInd in rInds:
    ax3[0].plot(z[rInd, :], Er[rInd, :], label='r = {:.0f} mm'.format(r[rInd, 0]*1e3))
    ax3[1].plot(z[rInd, :], Ez[rInd, :], label='r = {:.0f} mm'.format(r[rInd, 0]*1e3))
ax3[0].set_ylim([-Emax, Emax])
ax3[0].set_xlabel('z [m]')
ax3[0].set_ylabel('Er [V/m]')
ax3[0].legend()
ax3[0].grid()
ax3[1].set_ylim([-Emax, Emax])
ax3[1].set_xlabel('z [m]')
ax3[1].set_ylabel('Ez [V/m]')
ax3[1].grid()
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
    ax1[0].set_ylim([0, 0.55])
    plt.show(block=False)

input("Press Enter to continue...")
