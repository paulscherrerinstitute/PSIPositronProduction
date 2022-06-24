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
BUNCH_DOWNSAMPLING = 10

TARGET_Z = +41.   # [mm]
#
# RF_TYPE = 'field_map_PSI_Sband'
# TARGET_Z = +20.   # [mm]
# SOL_HOMOG_BZ = 1.5   # [T]

# TODO
# target.Ne = 1e4;	% no. of e- simulated
# target.Np = 4.37e10;	% no. of e+ required

AMD_FIELDMAP = 'RFTrack/YongkeTool_V2/field/field_map_HTS_5coils_Apr2022.dat'
R_APERTURE = 20e-3   # [m]
AMD_L_HALF = 96.5   # [mm]

TRACK_AFTER_AMD = True

RF_FIELDMAP_DIM = '3D'
RF_N_STRUCTURES = 4   # 11
RF_N_CELLS = 30
RF_N_PERIODS_PER_STRUCTURE = np.round(RF_N_CELLS/rfttools.RF_CLIC_CELLS_PER_PERIOD)
RF_L_STRUCTURE = RF_N_PERIODS_PER_STRUCTURE * bd.C/rfttools.RF_CLIC_FREQ  # [m]
RF_L_FLANGE = 0.05461   # [m]
RF_L_MECH_MARGIN = 0.1   # [m]
RF_T0 = 17.5 + 287.2 - 00./360.*150.   # [mm/c]
RF_PHASE_1ST = 171.   # [deg]
RF_PHASE_AFTER_1ST = 171.   # [deg]
RF_GRADIENT_1ST = 17.5e6   # [V/m]
RF_GRADIENT_AFTER_1ST = 21.0e6   # [V/m]
INITIAL_L = 0.   # [m]
RF_SEPARATION = 0.2   # [m]

# SOLENOID_TYPE = 'HomogeneousChannel'
SOL_HOMOG_BZ = 0.5   # [T]
#
SOLENOID_TYPE = 'Analytical'
SOL_R_IN_COIL = 0.115   # [m]
SOL_R_OUT_COIL = 0.250   # [m]
SOL_L = RF_L_STRUCTURE - RF_L_FLANGE - rfttools.RF_CLIC_L_CELL - RF_L_MECH_MARGIN
SOL_J = 3.15e6   # [A/m2]

FINAL_L = 0.   # [m]
T_MAX = 7500.   # [mm/c]

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
vol.set_aperture(R_APERTURE, R_APERTURE, 'circular')

amdFieldmap = bd.load_octave_matrices(AMD_FIELDMAP)
amdDz = amdFieldmap['Z'][0, 1] - amdFieldmap['Z'][0, 0]   # [mm]
amdDr = amdFieldmap['R'][1, 0]-amdFieldmap['R'][0, 0]   # [mm]
amdZ = amdFieldmap['Z'][0, :]   # [mm]
amdBzOnAxis = amdFieldmap['Bz'][0, :]

# Get effective field and length used in tracking,
# starts from target exit, stops at constant solenoid field value
indZTargetExit = (np.abs(amdZ - TARGET_Z)).argmin()
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
amdBrEff = amdFieldmap['Br'][:, indsEff]
amdBzEff = amdFieldmap['Bz'][:, indsEff]
amdZEff = amdZ[indsEff]
amdExitZInVolume = TARGET_EXIT_Z_IN_VOLUME + amdFieldLength*1e-3   # [m]
print('amdFieldLength = {:f} mm.'.format(amdFieldLength))

amd = rft.Static_Magnetic_FieldMap_2d(amdBrEff.T, amdBzEff.T, amdDr*1e-3, amdDz*1e-3)
if SOLENOID_TYPE == 'HomogeneousChannel':
    # TODO: Why is set_length() necessary?
    amd.set_length(amdFieldLength*1e-3)
    # TODO: Parametrize number of steps, i.e. following factor 2
    amd.set_nsteps(int(amdFieldLength*2.))
elif SOLENOID_TYPE == 'Analytical':
    amd.set_length(amdFieldLengthAnalytical*1e-3)
    amd.set_nsteps(int(amdFieldLengthAnalytical*2.))
# amd.set_aperture(R_APERTURE, R_APERTURE, 'circular')
amd.set_odeint_algorithm('rkf45')
vol.add(amd, 0, 0, TARGET_EXIT_Z_IN_VOLUME, 'entrance')
zFinalInVolume = amdExitZInVolume

if TRACK_AFTER_AMD:
    rf = rfttools.rf_clic_single_period(RF_FIELDMAP_DIM)
    lat = rft.Lattice()
    if INITIAL_L > 0:
        initialRfGap = rft.Drift(INITIAL_L)
        if SOLENOID_TYPE == 'HomogeneousChannel':
            initialRfGap.set_static_Bfield(0, 0, SOL_HOMOG_BZ)
        zFinalInVolume += INITIAL_L
        lat.append(initialRfGap)
    tRf = RF_T0
    rfGap = rft.Drift(RF_SEPARATION)
    if SOLENOID_TYPE == 'HomogeneousChannel':
        rf.set_static_Bfield(0, 0, SOL_HOMOG_BZ)
        rfGap.set_static_Bfield(0, 0, SOL_HOMOG_BZ)
    for structInd in np.arange(RF_N_STRUCTURES):
        lat.append(rfGap)
        # TODO: Get value form RF-Track object
        zFinalInVolume += RF_SEPARATION
        tRf += RF_SEPARATION * 1e3   # [mm/c]
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
        zFinalInVolume += RF_L_STRUCTURE / 2.56
    vol.add(lat, 0, 0, amdExitZInVolume, 'entrance')

if FINAL_L > 0:
    finalDrift = rft.Drift(FINAL_L)
    if SOLENOID_TYPE == 'HomogeneousChannel':
        finalDrift.set_static_Bfield(0, 0, SOL_HOMOG_BZ)
    vol.add(finalDrift, 0, 0, zFinalInVolume, 'entrance')
    zFinalInVolume += FINAL_L

vol.set_s0(0)
vol.set_s1(zFinalInVolume)

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
T0.t_max_mm = T_MAX
T0.tt_dt_mm = 1.   # [mm/c], track the emittance every tt_dt_mm (time)
T0.wp_dt_mm = 0.5e3   # [mm/c], save the distr. on disk every wp_dt_mm (time)
# Start tracking at s0
T0.backtrack_at_entrance = False
T0.odeint_algorithm = 'rk2'   # Options: 'rk2', 'rkf45', 'rk8pd'
T0.odeint_epsabs = 1e-5
T0.open_boundaries = False
T0.verbosity = 1   # 0 (default), 1 or 2

print('Tracking ends at s1 = {:f} m or at t_max = {:f} mm/c.'.format(zFinalInVolume, T0.t_max_mm))

B1_6dT = vol.track(B0_6dT, T0)
M1_6dT = B1_6dT.get_phase_space()
B1_6d = vol.get_bunch_at_s1()
M1_6d = B1_6d.get_phase_space("%x %xp %y %yp %t %Pc")
# A_AMD_LOSS = B_AMD_6dT.get_lost_particles();

# TODO: Qbunch
beamOut6dT, _ = bd.convert_rftrack_to_standard_df(
    rftrackDf=M1_6dT, rftrackDfFormat='rftrack_Px_S', pdgId=-11,
    outFwfPath='./DistrOut_AMD_6dT'
)
beamOut6d, _ = bd.convert_rftrack_to_standard_df(
    rftrackDf=M1_6d, rftrackDfFormat='rftrack_xp_t', pdgId=-11,
    outFwfPath='./DistrOut_AMD_6d'
)

fig1, ax1 = plt.subplots(5, 1)
rfttools.plot_transport(ax1, vol, B0_6dT, B1_6dT)
plt.show(block=False)
input("Press Enter to continue...")

# if (0)
#   V.set_s1(zv.HTS_exit*1e-3); %% [m]
#   B_AMD_6dT = V.track(B_target_6dT,T);
#   B_AMD_6d = V.get_bunch_at_s1();
#   A_HTS  = B_AMD_6d.get_phase_space("%x %xp %y %yp %t %Pc");
# endif

# % Yield

# amd.np_targ = rows(A_target);
# amd.np_amd  = rows(A_AMD);
# amd.efficiency = 1.0 * amd.np_amd / amd.np_targ;
# amd.yield      = 1.0 * amd.np_amd / target.Ne;

# printf(
#     "INFO_AMD:: AMD e+ collection efficiency: %.0f%% \n", amd.efficiency*100
# );
# printf("INFO_AMD:: AMD e+ yield: %.2f \n", amd.yield);

# MR = hypot(A_AMD(:,1), A_AMD(:,3)) < rf.R1i;
# amd.np_amd_r_acc = rows(A_AMD(MR,:));
# amd.efficiency_r_acc   = 1.0 * amd.np_amd_r_acc / amd.np_targ;
# amd.yield_r_acc = 1.0 * amd.np_amd_r_acc / target.Ne;
# printf("INFO_AMD:: AMD e+ collection efficiency (within RF R acceptance):
# %.0f%% \n", amd.efficiency_r_acc*100);
# printf("INFO_AMD:: AMD e+ yield (within RF R acceptance): %.2f \n", amd.yield_r_acc);

# amd = rmfield(amd,'field');


# outfname = ['amd_output/HTS_5coils_' linac_type '.dat'];
# %save('-text',outfname,'A_AMD','A_AMD_LOSS','A_HTS','amd','target');
# save('-text',outfname,'A_AMD');
