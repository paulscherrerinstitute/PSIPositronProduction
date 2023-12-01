import RF_Track as rft
import OctavePythonInterface as opi
import RFTrackTools as rfttools
import BeamDynamics as bd
import SimulationData as sd
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# INPUT Conventional FODO from 780 MeV #############################################################
# TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = 0.  # [m]
#   initialization just to enable output of beamline setup
# #
# TRACK_ONLY_REF_PART = False
# BUNCH_FILEPATH = 'Data/RFTrack/CaptureLinac/' + \
#     'PositronLinac_15RFStruct_LBandLargeR_2D_SolenoidsType1and2/' + \
#     'DistrOut_After2ndTracking_6d.sdf_txt'
# RFTRACK_FORMAT = 'rftrack_xp_t'
# BUNCH_PDGID = -11
# PARTICLE_CHARGE = +1   # [e], +1 = positrons, -1 = electrons
# PARTICLE_MASS = rft.electronmass   # [MeV/c/c]
# BUNCH_Z = 0.
# BUNCH_DOWNSAMPLING = 20
# FILTER_SPECS_MAIN_BUNCH = 'MainBunch'
# #
# TRACK_AFTER_MATCHING_1 = True
# #
# # RF_FIELDMAP = 'RunningSimulations/RFTrack/YongkeTool_V3/field/field_map_LargeR_Lband.dat'
# # RF_FIELDMAP_DIM = '1D'
# # RF_MOMENTUM_GAIN_PER_STRUCTURE = 59.35  # [MeV/c]
# # or
# RF_FIELDMAP = 'Data/Fieldmaps/pLinacF3_full44cells_YZplane_dy2mm_dz0p1L.dat'
# RF_FIELDMAP_DIM = '2D'
# RF_MOMENTUM_GAIN_PER_STRUCTURE = 59.35  # [MeV/c]
# #
# RF_N_PERIODS_PER_STRUCTURE = None
# RF_FIELDMAP_GRAD = 20e6  # [V/m]
# RF_N_STRUCTURES = 13  # 13 (1.54 GeV) + 22 (2.84 GeV)
# RF_L_STRUCTURE = 3.240  # [m]
# RF_R_APERTURE = 30e-3  # [m]
# RF_FREQ = 2e9  # [Hz]
# RF_PHASES = (0., )  # [deg]
# RF_SET_GRADIENTS = (20e6, 20e6, 20e6, 20e6, 20e6)  # [V/m]
# RF_SAMPLING_STEPS_PER_PERIOD = 360. / 5.
# #
# INITIAL_L = 0.  # [m]
# #
# MATCHING_1 = {
#     'Type': 'QuadMatching1',
#     'QuadLength': 0.500000,  # [m]
#     'QuadStrengths': [-0.193742, 0.560282, -0.600896, 0.501008, -0.442391],  # [1/m2]
#     'BrhoRef': 780. / PARTICLE_CHARGE,  # [MV/c]
#     'DriftLength': 3.000000,  # [m]
#     'RadialAperture': 45e-3,  # [m]
# }
# MATCHING_1_STRENGTH_TUNING = 0.9
# #
# FODO_1 = {
#     'Type': 'QuadFodo1',
#     'QuadLength': 0.350000,  # [m]
#     'QuadStrength': 0.824542,  # [1/m2]
#     'BrhoRef': 780. / PARTICLE_CHARGE,  # [MV/c]
#     'DriftLength': 4.051648,  # [m]
# }
# FODO_1_STRENGTH_TUNING = 0.85
# TODO: _TUNING paraameters far from 1 because setup does not exactly correspond to that determined
#  with Elegant. This was an oversight which could be easily ocrrected.
# # Best trial with 1D RF field map:
# # MATCHING_1_STRENGTH_TUNING = 0.75, FODO_1_STRENGTH_TUNING = 0.75
# # --> Capture efficiency = 0.925
# #
# FINAL_L = 1.  # [m]
####################################################################################################

# INPUT Conventional FODO from 735 MeV, with Chicane ###############################################
TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = 0.  # [m]
#   initialization just to enable output of beamline setup
LAT_R_APERTURE = 30e-3  # [m]
LAT_SAMPLING_STEPS_PER_MM = None  # [1/mm]
#
TRACK_ONLY_REF_PART = False
BUNCH_FILEPATH = 'RFTrackOutput/Baseline/' + \
    'PositronLinac_Chicane2m_After5RFStructs_14RFStructs_Positrons/' + \
    'DistrOut_After2ndTracking_6d.sdf_txt'
RFTRACK_FORMAT = 'rftrack_xp_t'
BUNCH_PDGID = -11
PARTICLE_CHARGE = +1   # [e], +1 = positrons, -1 = electrons
PARTICLE_MASS = rft.electronmass   # [MeV/c/c]
BUNCH_Z = 0.
BUNCH_DOWNSAMPLING = 1
FILTER_SPECS_MAIN_BUNCH = 'MainBunch'
#
TRACK_AFTER_MATCHING_1 = True
#
# RF_FIELDMAP = 'RunningSimulations/RFTrack/YongkeTool_V3/field/field_map_LargeR_Lband.dat'
# RF_FIELDMAP_DIM = '1D'
# RF_MOMENTUM_GAIN_PER_STRUCTURE = 59.35  # [MeV/c]
# or
# RF_FIELDMAP = 'Data/Fieldmaps/pLinacF3_full44cells_YZplane_dy2mm_dz0p1L.dat'
RF_FIELDMAP = 'Data/Fieldmaps/pLinacF3_full44cells_m2_YZplane_dy1mm_dz0p025L.dat'
RF_FIELDMAP_DIM = '2D'
RF_SMOOTH = 0
RF_MOMENTUM_GAIN_PER_STRUCTURE = 58.42  # [MeV/c]
#
RF_N_PERIODS_PER_STRUCTURE = None
RF_FIELDMAP_GRAD = 20e6  # [V/m]
RF_N_STRUCTURES = 14  # 14 (1.54 GeV) + 22 (2.84 GeV)
RF_L_STRUCTURE = 3.240  # [m]
RF_R_APERTURE = 30e-3  # [m]
RF_FREQ = 2e9  # [Hz]
RF_PHASES = (0., )  # [deg]
RF_SET_GRADIENTS = (20e6, 20e6, 20e6, 20e6, 20e6)  # [V/m]
RF_SAMPLING_STEPS_PER_PERIOD = 360. / 1.
#
INITIAL_L = 0.  # [m]
#
MATCHING_1 = {
    'Type': 'QuadMatching1',
    'QuadLength': 0.500000,  # [m]
    'QuadStrengths': [-0.153088, 0.564855, -0.691750, 0.428890, -0.441825],  # [1/m2]
    'BrhoRef': 735. / PARTICLE_CHARGE,  # [MV/c]
    'DriftLength': 2.000000,  # [m]
    'RadialAperture': 45e-3,  # [m]
}
MATCHING_1_STRENGTH_TUNING = 1.0
#
FODO_1 = {
    'Type': 'QuadFodo1',
    'QuadLength': 0.350000,  # [m]
    'QuadStrength': 0.8404469,  # [1/m2]
    'BrhoRef': 735. / PARTICLE_CHARGE,  # [MV/c]
    'DriftLength': 4.005464,  # [m]
}
FODO_1_STRENGTH_TUNING = 1.0
# Best trial:
#   MATCHING_1_STRENGTH_TUNING = , FODO_1_STRENGTH_TUNING =
#   --> Capture efficiency = 0.925
#
FINAL_L = 1.  # [m]
###################################################################################################


OUT_REL_PATH = 'RFTrackOutput/LatestSimPositronLinacSection3/'

M0 = bd.convert_standard_df_to_rftrack(
    sourceFilePath=BUNCH_FILEPATH, rftrackDfFormat=RFTRACK_FORMAT
)[0].to_numpy()
bunchPopulation = M0.shape[0]
# Definition of bunch population in RF-Track:
# total number of real particles in the bunch (not relevant here)
B0_6d = rft.Bunch6d(PARTICLE_MASS, bunchPopulation, PARTICLE_CHARGE, M0)

lat = rft.Lattice()
# TODO: Is the following not allowed as in Volume()?
# lat.set_aperture(LAT_R_APERTURE, LAT_R_APERTURE, 'circular')
beamlineSetup = pd.DataFrame(
    columns=['ElementType', 'zWrtTargetExit', 'MechanicalLength', 'Fieldmap'])

refPart1 = bd.get_json_entry(BUNCH_FILEPATH, FILTER_SPECS_MAIN_BUNCH, 'RefParticle1')
refPart1 = np.array([
    0., 0., 0., 0., refPart1['t']*bd.C*1e-6, refPart1['pz'], PARTICLE_MASS, PARTICLE_CHARGE, +1.])

zFinal = 0.  # [m]
if INITIAL_L > 0:
    initialDrift = rft.Drift(INITIAL_L)
    initialDrift.set_aperture(
        MATCHING_1['RadialAperture'], MATCHING_1['RadialAperture'], 'circular')
    lat.append(initialDrift)
    zFinal += initialDrift.get_length()

driftMatching1 = rft.Drift(MATCHING_1['DriftLength'])
driftMatching1.set_aperture(MATCHING_1['RadialAperture'], MATCHING_1['RadialAperture'], 'circular')
for quadStrength in MATCHING_1['QuadStrengths']:
    quadStrengthRftrack = quadStrength * MATCHING_1['QuadLength'] * MATCHING_1['BrhoRef'] \
            * MATCHING_1_STRENGTH_TUNING  # [MV/c/m]
    quadTmp = rft.Quadrupole(MATCHING_1['QuadLength'], quadStrengthRftrack)
    quadTmp.set_aperture(MATCHING_1['RadialAperture'], MATCHING_1['RadialAperture'], 'circular')
    if LAT_SAMPLING_STEPS_PER_MM is not None:
        quadTmp.set_nsteps(int(quadTmp.get_length() * 1e3 * LAT_SAMPLING_STEPS_PER_MM))
    lat.append(quadTmp)
    zFinal += quadTmp.get_length() / 2.
    beamlineSetup.loc[len(beamlineSetup.index)] = [
        'QuadMatching1', zFinal, MATCHING_1['QuadLength'],
        '{:.3f} MV/c/m (no fieldmap)'.format(quadTmp.get_strength())]
    zFinal += quadTmp.get_length() / 2.
    lat.append(driftMatching1)
    zFinal += driftMatching1.get_length()

quadPolarity = 1.
quadFodo1 = rft.Quadrupole(FODO_1['QuadLength'], 0)
if RF_R_APERTURE is not None:
    quadFodo1.set_aperture(RF_R_APERTURE, RF_R_APERTURE, 'circular')
if TRACK_AFTER_MATCHING_1:
    rfField = opi.load_octave_matrices(RF_FIELDMAP)
    # TODO: Refactorize this?
    if RF_N_PERIODS_PER_STRUCTURE is not None:
        rf = rfttools.rf_struct_from_single_period(
            rfField, RF_FIELDMAP_DIM, RF_N_PERIODS_PER_STRUCTURE, None, None, None,
            aperture=RF_R_APERTURE)
    else:
        rf = rfttools.rf_from_field_map(
            rfField, RF_FIELDMAP_DIM, None, None, None, smooth=RF_SMOOTH, aperture=RF_R_APERTURE)
    rf.set_nsteps(int(rf.get_length() / (bd.C/RF_FREQ) * RF_SAMPLING_STEPS_PER_PERIOD))
    BrhoRef = FODO_1['BrhoRef']
    quadRfGap = rft.Drift((FODO_1['DriftLength'] - rf.get_length()) / 2.)
    if RF_R_APERTURE is not None:
        quadRfGap.set_aperture(RF_R_APERTURE, RF_R_APERTURE, 'circular')
    for structInd in np.arange(RF_N_STRUCTURES):
        # Quad
        quadStrengthRftrack = FODO_1['QuadStrength'] * FODO_1['QuadLength'] * BrhoRef  # [MV/c/m]
        quadFodo1.set_strength(quadPolarity*quadStrengthRftrack*FODO_1_STRENGTH_TUNING)
        if LAT_SAMPLING_STEPS_PER_MM is not None:
            quadFodo1.set_nsteps(int(quadTmp.get_length() * 1e3 * LAT_SAMPLING_STEPS_PER_MM))
        lat.append(quadFodo1)
        zFinal += quadFodo1.get_length() / 2.
        beamlineSetup.loc[len(beamlineSetup.index)] = [
            'QuadFodo1', zFinal, FODO_1['QuadLength'],
            '{:.3f} MV/c/m (no fieldmap)'.format(quadFodo1.get_strength())
        ]
        zFinal += quadFodo1.get_length() / 2.
        quadPolarity *= -1
        # Drift
        lat.append(quadRfGap)
        zFinal += quadRfGap.get_length()
        # RF
        try:
            powerScalingFactor = (RF_SET_GRADIENTS[structInd] / RF_FIELDMAP_GRAD) ** 2.
        except IndexError:
            powerScalingFactor = (RF_SET_GRADIENTS[-1] / RF_FIELDMAP_GRAD) ** 2.
        rf.set_P_actual(powerScalingFactor)
        try:
            rfPhase = RF_PHASES[structInd]
        except IndexError:
            rfPhase = RF_PHASES[-1]
        rf.set_phid(rfPhase)
        lat.append(rf)
        zFinal += rf.get_length() / 2.
        beamlineSetup.loc[len(beamlineSetup.index)] = [
            'RF', zFinal, RF_L_STRUCTURE, os.path.basename(RF_FIELDMAP)]
        zFinal += rf.get_length() / 2.
        BrhoRef += RF_MOMENTUM_GAIN_PER_STRUCTURE
        # Drift
        lat.append(quadRfGap)
        zFinal += quadRfGap.get_length()
    print('Autophasing in lattice...')
    pzFinalAutophasing = lat.autophase(rft.Bunch6d(refPart1))
    print('pzFinalAutophasing = {:e} MeV/c'.format(pzFinalAutophasing))

if FINAL_L > 0:
    finalDrift = rft.Drift(FINAL_L)
    if RF_R_APERTURE is not None:
        finalDrift.set_aperture(RF_R_APERTURE, RF_R_APERTURE, 'circular')
    lat.append(finalDrift)
    zFinal += finalDrift.get_length()

beamlineSetup['zWrtAmdPeakField'] = \
    beamlineSetup['zWrtTargetExit'] + TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD
beamlineSetup[[
    'ElementType', 'zWrtTargetExit', 'zWrtAmdPeakField', 'MechanicalLength', 'Fieldmap'
]].to_csv(os.path.join(OUT_REL_PATH, 'BeamlineSetup.dat'), index=None)

# TODO: Verify tracking options with Andrea
# trackingOpts = rft.TrackingOptions()
lat.backtrack_at_entrance = False
lat.odeint_algorithm = 'rkf45'  # Options: 'rk2', 'rkf45', 'rk8pd'
lat.odeint_epsabs = 1e-5
lat.verbosity = 1  # 0 (default), 1 or 2
# TODO: Following parameters in Lattice()?
# trackingOpts.dt_mm = 0.2
# trackingOpts.tt_dt_mm = 1.  # [mm/c], track the emittance every tt_dt_mm (time)
# trackingOpts.wp_dt_mm = 0.5e3  # [mm/c], save the distr. on disk every wp_dt_mm (time)

if TRACK_ONLY_REF_PART:
    B1_6d = lat.track(rft.Bunch6dT(refPart1))
else:
    B1_6d = lat.track(B0_6d)
M1_6d = B1_6d.get_phase_space("%x %xp %y %yp %t %Pc")
bd.convert_rftrack_to_standard_df(
    rftrackDf=M1_6d, rftrackDfFormat='rftrack_xp_t', s=zFinal*1e3, pdgId=BUNCH_PDGID,
    outFwfPath=os.path.join(OUT_REL_PATH, 'DistrOut_6d'))

fig1, ax1 = plt.subplots(9, 1)
rfttools.save_plot_transport(ax1, lat, B0_6d, B1_6d, OUT_REL_PATH)
plt.show(block=False)

print("Continue debugging to close...")
print("... closed!")
