
  RF_Track;

  disp('INFO:: tracking new large-R TW RF cavities . . .');

  %% Load beam

  rf_input_file = 'Dat/AMDOutputPositrons_E6GeV_SpotSize0.5mm_EmittXY15um_ConvTarget5X0_HTS_0.5T.dat';

  load(rf_input_file); % TW / SW RF input in Octave plain text format

  %% t0 at RF entrance (auto phases)

  target.thickness = 17.5; % mm

  rf.t0 = target.thickness + amd.field_length; % mm/c
  
  printf("INFO_RF:: reference time at RF entrance: %.2f mm/c \n", rf.t0);
  
  Bunch.mass       = RF_Track.electronmass; % MeV/c/c
  Bunch.population = 0; %% number of particles per bunch (argument not used)
  Bunch.charge     = +1; %% units of e+
      
  B_AMD_6d = Bunch6d(Bunch.mass, Bunch.population, Bunch.charge, A_AMD);

  %% Load field

  rf.field = load('field/field_map_LargeR_Lband.dat').field;

  rf.Bc = amd.Bc; % Tesla
  rf.Ri = 30; % inner aperture radius in mm

  rf.P0_ref = [ 0 0 0 0 0 100 ]; % time will be reset later

  % Preliminary parameters
  if (amd.Bc == 0.50)
    rf.n_cav  = 5; % Nb of structures or cavities
    rf.Phs = [ -125.7 -127.8 -132.0 -102.9 -95.0 ]; % optimal
    % rf.Phs = [ -130 -130 -135 -75 -75 ]; % optimal & rounded for given reference particle
  elseif (amd.Bc == 1.50)
    rf.n_cav  = 5; % Nb of structures or cavities
    rf.Phs = [ -138.4 -142.5 -133.6 -79.0 -67.3 ]; % optimal
    % rf.Phs = [ -130 -130 -135 -75 -75 ]; % optimal & rounded
  else
    rf.n_cav  = 5; % Nb of structures or cavities
    rf.Phs = [ -125.7 -127.8 -132.0 -102.9 -95.0 ]; % same with 0.5 T
    % rf.Phs = [ -130 -130 -135 -75 -75 ]; % optimal & rounded
  endif

  rf.l_cav = 3.24e3; % structure full length in mm (including two ~120 mm drift tubes)
  rf.Grd = [ 20 * ones(1,rf.n_cav) ]; % dec. acc. ave. gradients in MV/m
  rf.fgap = 0; % front gap, distance to where 0.5 T starts in mm
  rf.cgap = 0; % cavity gap, distance between cavities in mm (already included in field map ~ 120 mm)
  rf.Gap = [ rf.fgap rf.cgap * ones(1,rf.n_cav-1) ]; % mm

  RF = RF_FieldMap_1d(  rf.field.E, ...  %% electric field, V/m
                        (rf.field.Z(2) - rf.field.Z(1))*1e-3, ... %% dz, m
                        range(rf.field.Z)*1e-3, ... %% l, m
			rf.field.frequency, ... %% f, Hz
			+1); % direction, usually +1, -1 in case of different sign conventions

  RF.set_static_Bfield(0,0,rf.Bc); % tesla
  RF.set_odeint_algorithm("rkf45");
  RF.set_nsteps(3000);
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

  % auto-extend gradients and phases
  if (rf.n_cav > length(rf.Grd))
    rf.Grd(length(rf.Grd)+1:rf.n_cav) = rf.Grd(end);
  endif
  if (rf.n_cav > length(rf.Phs))
    rf.Phs(length(rf.Phs)+1:rf.n_cav) = rf.Phs(end);
  endif

  printf("INFO:: RF Gradients = [ %s] MV/m \n",sprintf("%.1f ",rf.Grd));
  printf("INFO:: RF Phases = [ %s] degree \n",sprintf("%.0f ",rf.Phs));

  % reference particle
  rf.P0_ref = [ 0 0 0 0 rf.t0 100 ];

  code.optimise = 1; % save more / less info

  for i_rf = 1 : rf.n_cav
    RF_GAP.set_length(rf.Gap_ext(i_rf)*1e-3); %% [m]
    LAT.append(RF_GAP);
    RF.set_phid(rf.Phs(i_rf)); % degree
    RF.set_P_actual((rf.Grd(i_rf)/20.0)^2); % average gradient in MV/m
    LAT.append(RF);

    e_ap = LAT.autophase(Bunch6d(RF_Track.electronmass, 0, +1, rf.P0_ref));

    if(i_rf == 1)
      printf("INFO_RF:: Setting auto phase of Lattice() . . . \n");
      printf("          Reference particle: [ %s] \n", sprintf("%.2f ",rf.P0_ref));
    endif

    if(i_rf == rf.n_cav)
      printf("INFO_RF:: Max. energy gain of auto phase: %.1f MeV\n", (e_ap-rf.P0_ref(6)) );
    endif

    %% Track 1st RF
    if (!code.optimise && i_rf==1) %% frist structure
      printf("INFO_RF:: tracking 1st RF structure only . . . \n");
      tic
        B_RF_1 = LAT.track(B_AMD_6d);
      toc
      A_RF_1 = B_RF_1.get_phase_space("%x %xp %y %yp %t %Pc");
      rf.np_rf_1 = rows(A_RF_1);
      if (rows(A_RF_1)==0) rf.e_mean_rf_1 = 0;
      else 		   rf.e_mean_rf_1 = mean(A_RF_1(:,6)); endif
      printf("INFO_RF:: After 1st RF:\n");
      printf("          nPositrons: %i\n", rf.np_rf_1);
      printf("          Mean energy: %f MeV\n", rf.e_mean_rf_1);
    endif
    if (!code.optimise && i_rf==2) %% 2nd structure
      printf("INFO_RF:: tracking upto 2nd RF structure only . . . \n");
      tic
        B_RF_2 = LAT.track(B_AMD_6d);
      toc
      A_RF_2 = B_RF_2.get_phase_space("%x %xp %y %yp %t %Pc");
      rf.np_rf_2 = rows(A_RF_2);
      if (rows(A_RF_2)==0) rf.e_mean_rf_2 = 0;
      else 		   rf.e_mean_rf_2 = mean(A_RF_2(:,6)); endif
      printf("INFO_RF:: After 2nd RF:\n");
      printf("          nPositrons: %i\n", rf.np_rf_2);
      printf("          Mean energy: %f MeV\n", rf.e_mean_rf_2);
    endif
    if (!code.optimise && i_rf==3) %% 3rd structure
      printf("INFO_RF:: tracking upto 3rd RF structure only . . . \n");
      tic
        B_RF_3 = LAT.track(B_AMD_6d);
      toc
      A_RF_3 = B_RF_3.get_phase_space("%x %xp %y %yp %t %Pc");
      rf.np_rf_3 = rows(A_RF_3);
      if (rows(A_RF_3)==0) rf.e_mean_rf_3 = 0;
      else 		   rf.e_mean_rf_3 = mean(A_RF_3(:,6)); endif
      printf("INFO_RF:: After 3rd RF:\n");
      printf("          nPositrons: %i\n", rf.np_rf_3);
      printf("          Mean energy: %f MeV\n", rf.e_mean_rf_3);
    endif
  endfor

# Prepare Ez and Bz for plotting
zAxis = linspace(0, LAT.get_length(), 10000);   # [m]
Ez = [];
Bz = [];
for z = zAxis
    [E, B] = LAT.get_field(0, 0, z*1e3, 0);   # x,y,z,t (mm, mm/c)
    Ez(end+1) = E(3);
    Bz(end+1) = B(3);
endfor
emFields = [zAxis', Ez', Bz'];
outFile = fopen('Dat/LatestSim/EMFields.dat', 'w');
fprintf(outFile, 'z,Ez,Bz\n');
dlmwrite(outFile, emFields, '-append');
fclose(outFile);
figure(1)
plotyy(zAxis, Bz, zAxis, Ez)

  printf("INFO_RF:: tracking all RF structures . . . \n");
  tic
    B_RF = LAT.track(B_AMD_6d);  # Bunch6d(RF_Track.electronmass, 0, +1, rf.P0_ref)
  toc
  
  A_RF   =   B_RF.get_phase_space("%x %xp %y %yp %t %Pc");
  A_RF_TT = LAT.get_transport_table("%S %mean_E %sigma_E %sigma_t %N %emitt_x %emitt_y"); %% transport table along S

  ME = A_RF(:,6) > 50 & A_RF(:,6) < 350; % MeV

  rf.np_rf   = rows(A_RF);
  rf.np_rf_50_350 = rows(A_RF(ME,:));

  printf(" INFO_RF:: After all RFs:\n");
  printf("         nPositrons: %i\n", rf.np_rf);
  printf(" INFO_RF:: After all RFs with cut (50 < e < 350 MeV):\n");
  printf("         nPositrons: %i\n", rf.np_rf_50_350);

  if (1) % show mean energy

    if (rows(A_RF)==0) rf.e_mean_rf = 0;
    else 	       rf.e_mean_rf = mean(A_RF(:,6)); endif

    if (rows(A_RF(ME,6))==0) rf.e_mean_rf_truncated_50_350 = 0;
    else 	             rf.e_mean_rf_truncated_50_350 = mean(A_RF(ME,6)); endif

    printf("INFO_RF:: mean energy after all RFs: %f MeV\n", rf.e_mean_rf);
    printf("INFO_RF:: mean energy with cut (50 < e < 350 MeV) after all RFs: %f MeV\n", rf.e_mean_rf_truncated_50_350);
  endif

  %% Get total field for plotting
  if (!code.optimise && 1)
    %Z_pl  = linspace(V.get_s0()*1e3, V.get_s1()*1e3, 300); % mm
    Z_pl  = linspace(0, RF.get_length()*1e3, 1000); % mm
    Ez_pl = [];
    Bz_pl = [];
    for z = Z_pl
      [E_pl, B_pl] = RF.get_field(-1, 0, z, 0);
      Ez_pl = [ Ez_pl E_pl(3) ];
      Bz_pl = [ Bz_pl B_pl(3) ];
    end
    if (Ez_pl(end)==0)
      Ez_pl(end) = Ez_pl(end-1);
    endif
    if (Bz_pl(end)==0)
      Bz_pl(end) = rf.Bc;
    endif
    rf.Z_pl = Z_pl;
    rf.Ez_pl = Ez_pl;
    rf.Bz_pl = Bz_pl;
  endif

  rf = rmfield(rf,'field');

  outfname = rf_input_file;
  outfname = strrep(outfname,'AMDOutput','RF200MeVOutput');
  outfname = strrep(outfname,'.dat','_LargeRLband.dat');
  outfname = strrep(outfname,'Dat/','Dat/LatestSim/');
  if (code.optimise)
    save('-text',outfname,'A_RF','target','amd','rf');
  elseif (rf.n_cav >= 3)
    save('-text',outfname,'A_RF_1','A_RF_2','A_RF_3','A_RF','A_RF_TT','target','amd','rf');
  elseif (rf.n_cav == 2)
    save('-text',outfname,'A_RF_1','A_RF_2','A_RF','A_RF_TT','target','amd','rf');
  elseif (rf.n_cav == 1)
    save('-text',outfname,'A_RF_1','A_RF','A_RF_TT','target','amd','rf');
  endif

% TODO: Document difference between sigma_X (working with volume) and sigma_x (working with lattice)
strTT = '%mean_t %emitt_x %emitt_y %emitt_4d %sigma_x %sigma_y %mean_E';
TT = LAT.get_transport_table(strTT);
outFile = fopen('Dat/LatestSim/TransportTable.dat', 'w');
fprintf(outFile, [strrep(strrep(strTT,'%', ''),' ',',') '\n']);
dlmwrite(outFile, TT, '-append');
fclose(outFile);
tLims = [TT(1,1), TT(end,1)] / 1e3;
Elims = [0, TT(end,7)];
% Compute capture efficiency
Mlost = B_RF.get_lost_particles();   # Columns 1-6 like A_RF, z [mm] at which particle was lost, m [kg], Q [?] of particle type, Q of macro-particle [?]
% TODO: Check that index 7 is correct and that above comment is correct
if (~isempty(Mlost))
  Mlost = sortrows(Mlost, 7);
  sCapture = Mlost(:, 7);
  captureEff = (size(A_AMD,1) - (1:size(Mlost,1))') / size(A_AMD,1);
else
  sCapture = TT([1 end],1)/1e3;
  captureEff = [1.; 1.];
endif
captureEfficiency = [sCapture, captureEff];
outFile = fopen('Dat/LatestSim/CaptureEfficiency.dat', 'w');
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
