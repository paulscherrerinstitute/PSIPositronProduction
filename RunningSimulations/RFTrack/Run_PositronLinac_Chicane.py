import RF_Track as rft
import OctavePythonInterface as opi
import RFTrackTools as rfttools
import BeamDynamics as bd
import SimulationData as sd
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def append_drift(lat, zFinal, driftSpecs):
    if driftSpecs['L'] > 0:
        drift = rft.Drift(driftSpecs['L'])
        lat.append(drift)
        # drift.set_aperture(driftSpecs['RadialAperture'], driftSpecs['RadialAperture'], 'circular')
        zFinal += drift.get_length()
    return lat, zFinal


def append_dipole(lat, beamlineSetup, zFinal, dipoleSpecs):
    if dipoleSpecs['L'] > 0:
        dipole = rft.SBend(
            dipoleSpecs['L'], dipoleSpecs['BendingAngle'], dipoleSpecs['PoverQ'],
            dipoleSpecs['EdgeAngleIn'], dipoleSpecs['EdgeAngleOut'])
        # Drift with static By or RF field map with frequency 0 or set static_magnetic_field
        # (field only in volume, no yokes)
        lat.append(dipole)
        beamlineSetup.loc[len(beamlineSetup.index)] = [
            'Dipole1Chicane', zFinal, dipoleSpecs['L'],
            '{:.3f} MeV/c/C (no fieldmap)'.format(dipoleSpecs['PoverQ'])]
        # Set aperture
        zFinal += dipole.get_length()
    return lat, beamlineSetup, zFinal


# INPUT Chicane like SuperKEK-B ####################################################################
TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD = 0.035
#   [m], initialization just to enable output of beamline setup
#
TRACK_ONLY_REF_PART = False
BUNCH_FILEPATH = 'RFTrackOutput/LatestSimCaptureLinac/DistrOut_After1stTracking_6d.sdf_txt'
RFTRACK_FORMAT = 'rftrack_xp_t'
BUNCH_PDGID = -11
PARTICLE_CHARGE = +1   # [e], +1 = positrons, -1 = electrons
PARTICLE_MASS = rft.electronmass   # [MeV/c/c]
BUNCH_Z = 0.
BUNCH_DOWNSAMPLING = 10
FILTER_SPECS_MAIN_BUNCH = 'MainBunch'
#
CHICANE = {
    'AfterRFStructNo': 5,
    'Drift1': {
        'L': 0  # [m]
    },
    'Dipole1': {
        'L': 0.180,  # [m]
        'BendingAngle': 0.0526478,  # [rad]
        'PoverQ': 205.,  # [MeV/c/e]
        'EdgeAngleIn': 0,  # [rad]
        'EdgeAngleOut': 0,  # [rad]
    },
    'Drift2': {
        'L': 0.125  # [m]
    },
    'Dipole2': {
        'L': 0.180,  # [m]
        'BendingAngle': -0.0526478,  # [rad]
        'PoverQ': 205.,  # [MeV/c/e]
        'EdgeAngleIn': 0,  # [rad]
        'EdgeAngleOut': 0,  # [rad]
    },
    'Drift3': {
        'L': 0.125  # [m]
    },
    'Dipole3': {
        'L': 0.180,  # [m]
        'BendingAngle': -0.0526478,  # [rad]
        'PoverQ': 205.,  # [MeV/c/e]
        'EdgeAngleIn': 0,  # [rad]
        'EdgeAngleOut': 0,  # [rad]
    },
    'Drift4': {
        'L': 0.125  # [m]
    },
    'Dipole4': {
        'L': 0.180,  # [m]
        'BendingAngle': 0.0526478,  # [rad]
        'PoverQ': 205.,  # [MeV/c/e]
        'EdgeAngleIn': 0,  # [rad]
        'EdgeAngleOut': 0,  # [rad]
    },
    'Drift5': {
        'L': 0.125  # [m]
    },
    # 'TotLength': 1.5,  # [m]
}
####################################################################################################


OUT_REL_PATH = 'RFTrackOutput/LatestSimPositronLinacChicane/'

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
    columns=['ElementType', 'zWrtTargetExit', 'MechanicalLength', 'Fieldmap']
)

refPart1 = bd.get_json_entry(
    BUNCH_FILEPATH, FILTER_SPECS_MAIN_BUNCH, 'RefParticle1')
refPart1 = np.array([
    0., 0., 0., 0., refPart1['t']*bd.C*1e-6, refPart1['pz'], PARTICLE_MASS, PARTICLE_CHARGE, +1.
])

zFinal = 0.  # [m]

lat, zFinal = append_drift(lat, zFinal, CHICANE['Drift1'])
lat, beamlineSetup, zFinal = append_dipole(lat, beamlineSetup, zFinal,  CHICANE['Dipole1'])
lat, zFinal = append_drift(lat, zFinal, CHICANE['Drift2'])
lat, beamlineSetup, zFinal = append_dipole(lat, beamlineSetup, zFinal,  CHICANE['Dipole2'])
lat, zFinal = append_drift(lat, zFinal, CHICANE['Drift3'])
lat, beamlineSetup, zFinal = append_dipole(lat, beamlineSetup, zFinal,  CHICANE['Dipole3'])
lat, zFinal = append_drift(lat, zFinal, CHICANE['Drift4'])
lat, beamlineSetup, zFinal = append_dipole(lat, beamlineSetup, zFinal,  CHICANE['Dipole4'])
lat, zFinal = append_drift(lat, zFinal, CHICANE['Drift5'])

beamlineSetup['zWrtAmdPeakField'] = \
    beamlineSetup['zWrtTargetExit'] + TARGET_EXIT_Z_WRT_AMD_PEAK_FIELD
beamlineSetup[[
    'ElementType', 'zWrtTargetExit', 'zWrtAmdPeakField', 'MechanicalLength', 'Fieldmap'
]].to_csv(os.path.join(OUT_REL_PATH, 'BeamlineSetup.dat'), index=None)

# TODO: Verify tracking options with Andrea
trackingOpts = rft.TrackingOptions()
trackingOpts.tt_dt_mm = 1.   # [mm/c], track the emittance every tt_dt_mm (time)
# trackingOpts.wp_dt_mm = 0.5e3   # [mm/c], save the distr. on disk every wp_dt_mm (time)
trackingOpts.backtrack_at_entrance = False
trackingOpts.odeint_algorithm = 'rkf45'   # Options: 'rk2', 'rkf45', 'rk8pd'
trackingOpts.odeint_epsabs = 1e-5
trackingOpts.verbosity = 1  # 0 (default), 1 or 2

if TRACK_ONLY_REF_PART:
    B1_6d = lat.track(rft.Bunch6dT(refPart1))
else:
    B1_6d = lat.track(B0_6d)
M1_6d = B1_6d.get_phase_space("%x %xp %y %yp %t %Pc")
bd.convert_rftrack_to_standard_df(
    rftrackDf=M1_6d, rftrackDfFormat='rftrack_xp_t', s=zFinal*1e3, pdgId=BUNCH_PDGID,
    outFwfPath=os.path.join(OUT_REL_PATH, 'DistrOut_6d')
)

fig1, ax1 = plt.subplots(7, 1)
rfttools.save_plot_transport_in_lattice(ax1, lat, B0_6d, B1_6d, OUT_REL_PATH)
plt.show(block=False)

input("Press Enter to continue...")
