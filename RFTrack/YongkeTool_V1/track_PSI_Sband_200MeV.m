
  RF_Track;

  disp('INFO:: Tracking PSI S-band');

  %% Load beam

  rf_input_file = 'input/HTS_5coils_PSI_Sband.dat';

  load(rf_input_file); % TW / SW RF input in plain text format

  %% t0 at RF entrance

  rf.t0 = 17.5 + 179.4; % mm/c
  printf("INFO_RF:: reference time at RF entrance: %.2f mm/c \n", rf.t0);
  
  Bunch.mass       = RF_Track.electronmass; % MeV/c/c
  Bunch.population = 4.37E+10; %% number of particles per bunch
  Bunch.charge     = +1; %% units of e+

  B_AMD_6d = Bunch6d(Bunch.mass, Bunch.population, Bunch.charge, A_AMD);

  %% Load field

  rf.field = load('field/field_map_PSI_Sband.dat').Structure;

  rf.Bc = 1.5; % constant field in T
  rf.Ri = 20; % inner aperture radius in mm

  rf.n_cav  = 13; % Nb of structures or cavities
  rf.Phs = [ 220 220 215 215 210 210 220 270 285 270 240 235 235 ];

  rf.l_cav = 1.2e3; % structure full length in mm
  rf.Grd = [ 18 18 * ones(1,rf.n_cav-1) ]; % dec. acc. ave. gradients in MV/m
  rf.fgap = 0; % front gap, distance to where 0.5 T starts in mm
  rf.cgap = 150; % cavity gap, distance between cavities in mm
  rf.Gap = [ rf.fgap rf.cgap * ones(1,rf.n_cav-1) ]; % mm

  % 21 cells + couplers (1168 mm in total)
  % So each coupler is 59 mm long, with 16 mm margin.
  % One-Quarter field automatically extended to full range
  RF = RF_FieldMap(  rf.field.Ex, ...  %% electric field, V/m
                     rf.field.Ey, ...
                     rf.field.Ez, ...
                     rf.field.Bx, ...  %% magnetic field, T
                     rf.field.By, ...
                     rf.field.Bz, ...
                     rf.field.xa(1), ... %% m
                     rf.field.ya(1), ...
                     rf.field.xa(2) - rf.field.xa(1), ...
                     rf.field.ya(2) - rf.field.ya(1), ...
                     rf.field.za(2) - rf.field.za(1), ...
                     range(rf.field.za), ...
                     rf.field.frequency, rf.field.direction);

  RF.set_static_Bfield(0,0,rf.Bc); % tesla
  RF.set_cylindrical(false);
  RF.set_odeint_algorithm("rkf45");
  RF.set_nsteps(1000);
  RF.set_aperture(rf.Ri*1e-3, rf.Ri*1e-3, "circular"); % meter

  % regular cavity gap
  RF_GAP = Drift(rf.cgap*1e-3); % meter
  RF_GAP.set_static_Bfield(0,0,rf.Bc);
  RF_GAP.set_odeint_algorithm("rkf45");
  RF_GAP.set_nsteps(50);
  RF_GAP.set_aperture(0.1, 0.1, "circular"); % meter

  % structure margin included in gap
  rf.l_margin = (rf.l_cav - RF.get_length()*1e3)/2;
  rf.Gap_ext = rf.Gap;
  rf.Gap_ext(1) = rf.Gap(1) + rf.l_margin;
  rf.Gap_ext(2:end) = rf.Gap(2:end) + rf.l_margin*2;

  % lattice

  LAT = Lattice();

  for i_rf = 1 : rf.n_cav
    RF_GAP.set_length(rf.Gap_ext(i_rf)*1e-3); %% [m]
    LAT.append(RF_GAP);
    RF.set_phid(rf.Phs(i_rf)); % degree
    RF.set_P_actual((rf.Grd(i_rf)/18.0)**2); % average gradient in MV/m
    RF.set_t0( rf.t0 + sum(rf.Gap(1:i_rf)) + rf.l_cav*(i_rf-1) ); % mm/c
    LAT.append(RF);
  endfor

  printf("INFO_RF:: tracking all RF structures . . . \n");
  tic
    B_RF = LAT.track(B_AMD_6d);
  toc
  
  A_RF   =   B_RF.get_phase_space("%x %xp %y %yp %t %Pc");

  save('-text','output/HTS_5coils_PSI_Sband.dat','A_RF');
