
  RF_Track;

  p_type = 'e+';

  linac_type = 'CLIC_Lband';
  %linac_type = 'PSI_Sband';

  disp(['INFO:: Particle: ' p_type]);
  disp(['INFO:: Capture Linac: ' linac_type]);

  %%

  target.Ne = 1e4;	% no. of e- simulated
  target.Np = 4.37e10;	% no. of e+ required

  %%

  amd.fieldmap_file = 'field/field_map_HTS_5coils_Apr2022.dat';

  amd.Ri  = 20; % Effective inner aperture radius in mm

  amd.half_length = 96.5; % half length of HTS solenoid in mm

  if (strcmp(linac_type,'CLIC_Lband'))
    amd.Ztc = +41; % Target exit location w.r.t HTS upstream solenoid exit in mm
    amd.Bc  = 0.5; % Constant solenoid field for acceleration in T
  endif

  if (strcmp(linac_type,'PSI_Sband'))
    amd.Ztc = +20; % Target exit location w.r.t HTS upstream solenoid exit in mm
    amd.Bc  = 1.5; % Constant solenoid field for acceleration in T
  endif

  %%

  rf.R1i  = 20; % mm
  
  % Input
  
  amd_input_file = 'amd_input/E6GeV_SpotSize0.5mm_Target5X0.dat';
  
  try
    A_target = load(amd_input_file);
  catch
    printf("ERROR:: Empty AMD input file ! ! !\n");
    return;
  end_try_catch

  %

  %M = A_target(:,5) < 1000; % mm/c, remove very large time particles
  %A_target = A_target(M,:);
    
  Bunch.mass       = RF_Track.electronmass; % MeV/c/c
  Bunch.population = target.Np; %% number of particles per bunch

  if (strcmp(p_type,'e+'))
    Bunch.charge     = +1; %% units of e+
  endif

  if (strcmp(p_type,'e-'))
    Bunch.charge     = -1; %% units of e+
  endif
        
  B_target_6d = Bunch6d(Bunch.mass, Bunch.population, Bunch.charge, A_target);

  %%

  zv.target_exit = 0; % target exit position in Volume in mm, always being 0

  %% Load AMD field map

  amd.field = load(amd.fieldmap_file);
  
  amd.field.dz = amd.field.Z(1,2)-amd.field.Z(1,1); %% [mm]
  amd.field.dr = amd.field.R(2,1)-amd.field.R(1,1); %% [mm]

  amd.Z_axis  = amd.field.Z(1,:); % on-axis Z
  amd.Bz_axis = amd.field.Bz(1,:); % on-axis Bz

  % Get effective field and length used in tracking, starts from target exit, stops at constant solenoid field value

  amd.zfte = amd.Ztc; % target exit position in field map in mm

  izte = lookup(amd.Z_axis, amd.zfte);
  izbc = length(amd.Z_axis) - lookup(flip(amd.Bz_axis),amd.Bc) + 1;

  amd.field.Br_eff = amd.field.Br(:, izte:izbc);
  amd.field.Bz_eff = amd.field.Bz(:, izte:izbc);
  amd.Z_axis_eff  = amd.Z_axis(izte:izbc);
  amd.Bz_axis_eff = amd.Bz_axis(izte:izbc);

  amd.field_length = range(amd.Z_axis_eff); 

  %% AMD field

  HTS_FIELD = Static_Magnetic_FieldMap_2d(amd.field.Br_eff', amd.field.Bz_eff', amd.field.dr*1e-3, amd.field.dz*1e-3); %% [T, m]
  HTS_FIELD.set_length(amd.field_length*1e-3); %% [m]
  HTS_FIELD.set_nsteps(floor(amd.field_length*2));
  HTS_FIELD.set_odeint_algorithm("rkf45");

  %% HTS solenoid shape
  
  HTS_SHAPE = AdiabaticMatchingDevice( amd.half_length*2*1e-3, 0, 0); % [L, B0, Mu] in [m, T, 1/m]
  HTS_SHAPE.set_nsteps(200);
  HTS_SHAPE.set_odeint_algorithm("rk2");
  HTS_SHAPE.set_entrance_aperture(amd.Ri*1e-3); % m
  HTS_SHAPE.set_exit_aperture(amd.Ri*1e-3); % m
  HTS_SHAPE.set_static_Bfield(0,0,0); % T

  zv.HTS_exit = zv.target_exit - amd.Ztc + amd.half_length; % HTS exit position in Volume in mm

  %% Volume definition

  V = Volume();
    V.set_aperture(1, 1, 'circular'); % m
    V.add(HTS_SHAPE, 0, 0, zv.HTS_exit*1e-3, 'exit'); % m
    V.add(HTS_FIELD, 0, 0, zv.target_exit*1e-3, 'entrance'); % m
    V.set_s0(0); %% [m]
    V.set_s1(zv.target_exit*1e-3 + amd.field_length*1e-3); %% [m]

  %% Get total field for plotting
  if (1)
    Z_pl  = linspace(V.get_s0()*1e3, V.get_s1()*1e3, 300); % mm
    Bz_pl = [];
    for z = Z_pl
      [E_pl, B_pl] = V.get_field(0, 0, z, 0);
      Bz_pl = [ Bz_pl B_pl(3) ];
    end
    if (Bz_pl(end)==0)
      Bz_pl(end) = amd.Bc;
    endif
    Z_pl += zv.target_exit + amd.zfte; % in field map coordinate
    amd.Z_pl = Z_pl;
    amd.Bz_pl = Bz_pl;
  endif

  %% Track with Volume

  T = TrackingOptions();
    T.dt_mm = 0.1; %% [mm/c]
    T.t_max_mm = Inf;
    T.backtrack_at_entrance = false; %% start tracking at s0
    T.odeint_algorithm = "rkf45"; %% 'rk2', 'rkf45', 'rk8pd'
    T.odeint_epsabs = 1e-5;
    T.open_boundaries = false;

  disp('INFO:: tracking HTS with Volume . . .')
  tic;
    B_target_6dT = Bunch6dT(B_target_6d);
    B_AMD_6dT = V.track(B_target_6dT,T);
    B_AMD_6d = V.get_bunch_at_s1();
  toc;
  A_AMD  = B_AMD_6d.get_phase_space("%x %xp %y %yp %t %Pc");
  A_AMD_LOSS = B_AMD_6dT.get_lost_particles();

  if (0)
    V.set_s1(zv.HTS_exit*1e-3); %% [m]
    B_AMD_6dT = V.track(B_target_6dT,T);
    B_AMD_6d = V.get_bunch_at_s1();
    A_HTS  = B_AMD_6d.get_phase_space("%x %xp %y %yp %t %Pc");
  endif

  %% Yield

  amd.np_targ = rows(A_target);
  amd.np_amd  = rows(A_AMD);
  amd.efficiency = 1.0 * amd.np_amd / amd.np_targ;
  amd.yield      = 1.0 * amd.np_amd / target.Ne;

  printf("INFO_AMD:: AMD e+ collection efficiency: %.0f%% \n", amd.efficiency*100);
  printf("INFO_AMD:: AMD e+ yield: %.2f \n", amd.yield);

  MR = hypot(A_AMD(:,1), A_AMD(:,3)) < rf.R1i;
  amd.np_amd_r_acc = rows(A_AMD(MR,:));
  amd.efficiency_r_acc   = 1.0 * amd.np_amd_r_acc / amd.np_targ;
  amd.yield_r_acc = 1.0 * amd.np_amd_r_acc / target.Ne;
  printf("INFO_AMD:: AMD e+ collection efficiency (within RF R acceptance): %.0f%% \n", amd.efficiency_r_acc*100);
  printf("INFO_AMD:: AMD e+ yield (within RF R acceptance): %.2f \n", amd.yield_r_acc);

  amd = rmfield(amd,'field');


  outfname = ['amd_output/HTS_5coils_' linac_type '.dat'];
  %save('-text',outfname,'A_AMD','A_AMD_LOSS','A_HTS','amd','target');
  save('-text',outfname,'A_AMD');
