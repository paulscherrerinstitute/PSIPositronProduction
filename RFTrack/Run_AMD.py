import RF_Track as rft
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

RF_TYPE = 'field_map_CLIC_Lband'
TARGET_Z = +41.   # [mm]
BZ_SOLENOID = 0.5   # [T]
#
# RF_TYPE = 'field_map_PSI_Sband'
# TARGET_Z = +20.   # [mm]
# BZ_SOLENOID = 1.5   # [T]

# TODO
# target.Ne = 1e4;	% no. of e- simulated
# target.Np = 4.37e10;	% no. of e+ required

AMD_FIELDMAP = 'RFTrack/YongkeTool_V2/field/field_map_HTS_5coils_Apr2022.dat'
R_APERTURE = 20e-3   # [mm]
AMD_L_HALF = 96.5   # [mm]

# TODO
# rf.R1i  = 20; % mm

# TODO: Refactorize, same code in Run_Linac1_Section1_Simple.py
# A_target = load(BUNCH_FILEPATH)
distrMatNp = np.loadtxt(BUNCH_FILEPATH, skiprows=1)
beamIn, _ = bd.convert_rftrack_to_standard_df(
    rftrackDf=distrMatNp, rftrackDfFormat=RFTRACK_FORMAT,
    z=BUNCH_Z, pdgId=BUNCH_PDGID
)
# TODO:
# beamIn['Q'] = Q_DRIVE_BEAM / N_MACROPARTICLES_DRIVE_BEAM
# TODO:
# M = A_target(:,5) < 1000; % mm/c, remove very large time particles
# A_target = A_target(M,:);
# if FILTER_SPECS_SELECTOR is not None:
#     with open(sd.build_data_path(REL_PATH, 'filterSpecs.json'), 'r') as filterSpecsFile:
#         filterSpecsList = json.load(filterSpecsFile)
#     filterSpecs = filterSpecsList[FILTER_SPECS_SELECTOR]['filterSpecs']
#     refParticle = filterSpecsList[FILTER_SPECS_SELECTOR]['RefParticle1']
#     beamIn = bd.filter_distr(beamIn, filterSpecs)
M0 = bd.convert_standard_df_to_rftrack(standardDf=beamIn, rftrackDfFormat=RFTRACK_FORMAT)[0].to_numpy()
BUNCH_POPULATION = M0.shape[0]
B06d = rft.Bunch6d(PARTICLE_MASS, BUNCH_POPULATION, PARTICLE_CHARGE, M0)
B06dT = rft.Bunch6dT(B06d)


# zv.target_exit = 0; % target exit position in Volume in mm, always being 0
amdFieldmap = bd.load_octave_matrices(AMD_FIELDMAP)
amdDz = amdFieldmap['Z'][0,1] - amdFieldmap['Z'][0,0]   # [mm]
amdDr = amdFieldmap['R'][1,0]-amdFieldmap['R'][0,0]   # [mm]
amdZ  = amdFieldmap['Z'][0,:]
amdBzOnAxis = amdFieldmap['Bz'][0,:]

# % Get effective field and length used in tracking, starts from target exit, stops at constant solenoid field value
# amd.zfte = amd.Ztc; % target exit position in field map in mm
# izte = lookup(amd.Z_axis, amd.zfte);
# izbc = length(amd.Z_axis) - lookup(flip(amd.Bz_axis),amd.Bc) + 1;
# amd.field.Br_eff = amd.field.Br(:, izte:izbc);
# amd.field.Bz_eff = amd.field.Bz(:, izte:izbc);
# amd.Z_axis_eff  = amd.Z_axis(izte:izbc);
# amd.Bz_axis_eff = amd.Bz_axis(izte:izbc);
# amd.field_length = range(amd.Z_axis_eff); 

# % AMD field

# HTS_FIELD = Static_Magnetic_FieldMap_2d(amd.field.Br_eff', amd.field.Bz_eff', amd.field.dr*1e-3, amd.field.dz*1e-3); %% [T, m]
# HTS_FIELD.set_length(amd.field_length*1e-3); %% [m]
# HTS_FIELD.set_nsteps(floor(amd.field_length*2));
# HTS_FIELD.set_odeint_algorithm("rkf45");

# % HTS solenoid shape

# HTS_SHAPE = AdiabaticMatchingDevice( amd.half_length*2*1e-3, 0, 0); % [L, B0, Mu] in [m, T, 1/m]
# HTS_SHAPE.set_nsteps(200);
# HTS_SHAPE.set_odeint_algorithm("rk2");
# HTS_SHAPE.set_entrance_aperture(amd.Ri*1e-3); % m
# HTS_SHAPE.set_exit_aperture(amd.Ri*1e-3); % m
# HTS_SHAPE.set_static_Bfield(0,0,0); % T

# zv.HTS_exit = zv.target_exit - amd.Ztc + amd.half_length; % HTS exit position in Volume in mm

# % Volume definition

# V = Volume();
#   V.set_aperture(1, 1, 'circular'); % m
#   V.add(HTS_SHAPE, 0, 0, zv.HTS_exit*1e-3, 'exit'); % m
#   V.add(HTS_FIELD, 0, 0, zv.target_exit*1e-3, 'entrance'); % m
#   V.set_s0(0); %% [m]
#   V.set_s1(zv.target_exit*1e-3 + amd.field_length*1e-3); %% [m]

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

# T = TrackingOptions();
#   T.dt_mm = 0.1; %% [mm/c]
#   T.t_max_mm = Inf;
#   T.backtrack_at_entrance = false; %% start tracking at s0
#   T.odeint_algorithm = "rkf45"; %% 'rk2', 'rkf45', 'rk8pd'
#   T.odeint_epsabs = 1e-5;
#   T.open_boundaries = true;

# disp('INFO:: tracking HTS with Volume . . .')
# tic;
#   B_target_6dT = Bunch6dT(B_target_6d);
#   B_AMD_6dT = V.track(B_target_6dT,T);
#   B_AMD_6d = V.get_bunch_at_s1();
# toc;
# A_AMD  = B_AMD_6d.get_phase_space("%x %xp %y %yp %t %Pc");
# A_AMD_LOSS = B_AMD_6dT.get_lost_particles();

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

# printf("INFO_AMD:: AMD e+ collection efficiency: %.0f%% \n", amd.efficiency*100);
# printf("INFO_AMD:: AMD e+ yield: %.2f \n", amd.yield);

# MR = hypot(A_AMD(:,1), A_AMD(:,3)) < rf.R1i;
# amd.np_amd_r_acc = rows(A_AMD(MR,:));
# amd.efficiency_r_acc   = 1.0 * amd.np_amd_r_acc / amd.np_targ;
# amd.yield_r_acc = 1.0 * amd.np_amd_r_acc / target.Ne;
# printf("INFO_AMD:: AMD e+ collection efficiency (within RF R acceptance): %.0f%% \n", amd.efficiency_r_acc*100);
# printf("INFO_AMD:: AMD e+ yield (within RF R acceptance): %.2f \n", amd.yield_r_acc);

# amd = rmfield(amd,'field');


# outfname = ['amd_output/HTS_5coils_' linac_type '.dat'];
# %save('-text',outfname,'A_AMD','A_AMD_LOSS','A_HTS','amd','target');
# save('-text',outfname,'A_AMD');
