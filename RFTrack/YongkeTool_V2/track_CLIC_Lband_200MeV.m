
  RF_Track;

  disp('INFO:: tracking CLIC-like TW RF cavities . . .');

  %% Load beam

  % rf_input_file = 'rf_input/HTS_5coils_CLIC_Lband.dat';
  rf_input_file = 'amd_output/HTS_5coils_CLIC_Lband.dat';

  load(rf_input_file); % TW / SW RF input in plain text format

  %% t0 at RF entrance

  rf.t0 = 17.5 + 287.2; % mm/c
  printf("INFO_RF:: reference time at RF entrance: %.2f mm/c \n", rf.t0);
  
  Bunch.mass       = RF_Track.electronmass; % MeV/c/c
  Bunch.population = 4.37E+10; %% number of particles per bunch
  Bunch.charge     = +1; %% units of e+
      
  B_AMD_6d = Bunch6d(Bunch.mass, Bunch.population, Bunch.charge, A_AMD);

  %% Load field

  rf.field = load('field/field_map_CLIC_Lband.dat').M;

  rf.Bc = 0.5; % Tesla
  rf.Ri = 20; % inner aperture radius in mm
  rf.freq = 1.9986163867e+09; % in Hz
  rf.n_cell = 30; % Nb of cells per structure
  rf.n_cav  = 11; % Nb of structures or cavities (default 11)
  rf.n_wave = round(rf.n_cell/3); % Nb of waves per structure
  rf.Grd = [ 17.5 21.0 ]; % dec. acc. ave. gradients in MV/m
  rf.Phs = [ 171 171 ]; % dec. acc. phases in deg [default]
  % rf.Phs = [ 171+180 171+180 ]; % dec. acc. phases in deg [for test only]
  rf.fgap = 0; % front gap, distance to where 0.5 T starts in mm
  rf.cgap = 200; % cavity gap, distance between cavities in mm
  rf.Gap = [ rf.fgap ones(1,rf.n_cav-1)*rf.cgap ]; % al gaps in mm

  % 3 cells field (1/10 of full structure)
  RF = RF_FieldMap(  rf.field.Ex, ...  %% electric field, V/m
                     rf.field.Ey, ...
                     rf.field.Ez, ...
                     rf.field.Bx, ...  %% magnetic field, T
                     rf.field.By, ...
                     rf.field.Bz, ...
                     rf.field.ra(1), ...
                     rf.field.ta(1), ...
                     rf.field.ra(2) - rf.field.ra(1), ...
                     rf.field.ta(2) - rf.field.ta(1), ...
                     rf.field.za(2) - rf.field.za(1), ...
                     rf.field.za(end), ...
                     rf.freq, +1);

  RF.set_static_Bfield(0,0,rf.Bc); % tesla
  RF.set_cylindrical(true);
  RF.set_odeint_algorithm("rkf45");
  RF.set_nsteps(50);
  RF.set_aperture(20e-3, 20e-3, "circular"); % meter

  % regular cavity gap
  RF_GAP = Drift(rf.cgap*1e-3); % meter
  RF_GAP.set_static_Bfield(0,0,rf.Bc);
  RF_GAP.set_odeint_algorithm("rkf45");
  RF_GAP.set_nsteps(50);
  RF_GAP.set_aperture(0.1, 0.1, "circular"); % meter

  % lattice

  LAT = Lattice();

  for i_rf = 1 : rf.n_cav
    RF_GAP.set_length(rf.Gap(i_rf)*1e-3); %% [m]
    LAT.append(RF_GAP);
    if i_rf == 1
      % decelerating
      RF.set_phid(rf.Phs(1)); % degree
      RF.set_P_actual((rf.Grd(1)/11.23)^2); % average gradient in MV/m
    else
      % accelerating
      RF.set_phid(rf.Phs(2)); % degree
      RF.set_P_actual((rf.Grd(2)/11.23)^2); % average gradient in MV/m
    endif
    RF.set_t0( rf.t0 + sum(rf.Gap(1:i_rf)) ); % mm/c
    fprintf('Setting rf.t0 = %.1f mm/c\n', RF.get_t0())
    for i_wave = 1:rf.n_wave
      r_rf = rf.Ri;
      RF.set_aperture(r_rf*1e-3, r_rf*1e-3, "circular"); % meter
      LAT.append(RF);
    endfor
  endfor

# Prepare Ez and Bz for plotting
zAxis = linspace(0, LAT.get_length(), 1000);   # [m]
Ez = [];
Bz = [];
for z = zAxis
    [E, B] = LAT.get_field(0, 0, z*1e3, 0);   # x,y,z,t (mm, mm/c)
    Ez(end+1) = E(3);
    Bz(end+1) = B(3);
endfor
emFields = [zAxis', Ez', Bz'];
% save('-text', 'rf_output/LatestSim/EMFields.dat', 'emFields');
outFile = fopen('rf_output/LatestSim/EMFields.dat', 'w');
fprintf(outFile, 'z,Ez,Bz\n');
dlmwrite(outFile, emFields, '-append');
fclose(outFile);
figure(1)
%plot(zAxis, Bz)
plotyy(zAxis, Bz, zAxis, Ez)

  printf("INFO_RF:: tracking all RF structures . . . \n");
  T = TrackingOptions();
  T.tt_dt_mm = 1.;   % [mm], track the emittance every tt_dt_mm (time)
  T.verbosity = 1;
  tic
    B_RF = LAT.track(B_AMD_6d);
  toc
  
  A_RF   =   B_RF.get_phase_space("%x %xp %y %yp %t %Pc");

  save('-text','rf_output/LatestSim/HTS_5coils_CLIC_Lband.dat','A_RF');
% TODO: Document difference between sigma_X (working with volume) and sigma_x (working with lattice)
strTT = '%mean_t %emitt_x %emitt_y %emitt_4d %sigma_x %sigma_y %mean_E';
TT = LAT.get_transport_table(strTT);
% save('-text', 'rf_output/LatestSim/TransportTable.dat', 'TT');
outFile = fopen('rf_output/LatestSim/TransportTable.dat', 'w');
fprintf(outFile, [strrep(strrep(strTT,'%', ''),' ',',') '\n']);
dlmwrite(outFile, TT, '-append');
fclose(outFile);
tLims = [TT(1,1), TT(end,1)] / 1e3;
Elims = [0, TT(end,7)];
% Compute capture efficiency
Mlost = B_RF.get_lost_particles();   # Columns 1-6 like A_RF, z [mm] at which particle was lost, m [kg], Q [?] of particle type, Q of macro-particle [?]
% TODO: Check that index 7 is correct and that above comment is correct
Mlost = sortrows(Mlost, 7);
sCapture = Mlost(:, 7);
captureEff = (size(A_AMD,1) - (1:size(Mlost,1))') / size(A_AMD,1);
captureEfficiency = [sCapture, captureEff];
% save('-text', 'rf_output/LatestSim/CaptureEfficiency.dat', 'captureEfficiency');
outFile = fopen('rf_output/LatestSim/CaptureEfficiency.dat', 'w');
fprintf(outFile, 's,CaptureEfficiency\n');
dlmwrite(outFile, captureEfficiency, '-append');
fclose(outFile);
figure(2)
subplot(5, 1, 1)
plot(TT(:,1)/1e3, TT(:,7))
xlim(tLims)
ylim(Elims)
xlabel('s [m]')
ylabel('Beam energy [MeV]')
grid()
subplot(4, 1, 2)
plot(sCapture/1e3, captureEff)
xlim(tLims)
ylim([0, 1])
xlabel('s [m]')
ylabel('Capture efficiency')
grid()
subplot(4, 1, 3)
plot(TT(:,1)/1e3, TT(:,2))
hold on;
plot(TT(:,1)/1e3, TT(:,3))
plot(TT(:,1)/1e3, TT(:,4))
xlim(tLims)
xlabel('s [m]')
ylabel('emitt [pi mm mrad]')
legend('emitt x', 'emitt y', 'emitt 4d')
grid()
subplot(4, 1, 4)
plot(TT(:,1)/1e3, TT(:,5))
hold on
plot(TT(:,1)/1e3, TT(:,6))
xlim(tLims)
xlabel('s [m]')
ylabel('sigma [mm]')
legend('sigma x', 'sigma y')
grid()
waitforbuttonpress()
