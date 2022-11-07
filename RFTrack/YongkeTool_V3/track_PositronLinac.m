
  use_analytic_sol = 0;

  RF_Track;

  disp('INFO:: tracking e+ linac (analytic, longitudinal, 1.54 GeV) . . .');

  %% Load beam

  %%pl_input_file = 'Dat/RF200MeVOutputPositrons_E6GeV_SpotSize0.5mm_EmittXY15um_ConvTarget5X0_HTS_0.5T_CLICLband.dat';
  %%pl_input_file = 'Dat/RF200MeVOutputPositrons_E6GeV_SpotSize0.5mm_EmittXY15um_ConvTarget5X0_HTS_0.5T_PSISband.dat';
  % pl_input_file = 'Dat/RF200MeVOutputPositrons_E6GeV_SpotSize0.5mm_EmittXY15um_ConvTarget5X0_HTS_0.5T_LargeRLband.dat';
  % outfname = 'Dat/PositronLinacOutputPositrons_E6GeV_SpotSize0.5mm_EmittXY15um_ConvTarget5X0_HTS_0.5T_LargeRLband.dat'
  % pl_input_file = '../../Data/RFTrackResults/YongkeTool_V3/CaptureLinacUpTo200MeV_LBandLargeR/RF200MeVOutputPositrons_E6GeV_SpotSize0.5mm_EmittXY15um_ConvTarget5X0_HTS_0.5T_LargeRLband.dat';
  % outfname = '../../Data/RFTrackResults/YongkeTool_V3/CaptureLinacUpTo200MeV_LBandLargeR/PositronLinacOutputPositrons_E6GeV_SpotSize0.5mm_EmittXY15um_ConvTarget5X0_HTS_0.5T_LargeRLband.dat';
  % pl_input_file = '../../Data/RFTrackResults/CaptureLinac/CaptureLinacUpTo200MeV_LBandLargeR/DistrOut_After1stTracking_6d.dat';
  % outfname = '../../Data/RFTrackResults/CaptureLinac/CaptureLinacUpTo200MeV_LBandLargeR/DistrOut_PositronLinac.dat'
  % pl_input_file = '../../Data/RFTrackResults/CaptureLinac/CaptureLinacUpTo200MeV_LBandLargeR_RealSolenoids_Type1and2/DistrOut_After1stTracking_6d.dat';
  % outfname = '../../Data/RFTrackResults/CaptureLinac/CaptureLinacUpTo200MeV_LBandLargeR_RealSolenoids_Type1and2/DistrOut_PositronLinac.dat'
  pl_input_file = '../../Data/RFTrackResults/CaptureLinac/CaptureLinacUpTo200MeV_LBandLargeR_RealSolenoids_Type1and2_TargetAt35mm/DistrOut_After1stTracking_6d.dat';
  outfname = '../../Data/RFTrackResults/CaptureLinac/CaptureLinacUpTo200MeV_LBandLargeR_RealSolenoids_Type1and2_TargetAt35mm/DistrOut_PositronLinac.dat'

  if(use_analytic_sol)
  %  pl_input_file = 'Dat/AMDAndRF200MeVOutputPositrons_E6GeV_SpotSize0.5mm_EmittXY15um_ConvTarget5X0_LargeRLband_AnalyticSolenoid.dat';
  endif

  try
    % load(pl_input_file); % e+ linac input in plain text format
    A_RF = load(pl_input_file);
    printf(pl_input_file)
  catch
    printf("ERROR_PL:: Empty positron linac input file ! ! !\n");
    return;
  end_try_catch

  %%pl.freq = 2.856e9; % Hz
  pl.freq = 2e9; % Hz
  pl.e_target = 1.54; % GeV
  pl.e_accept = 3.8e-2; % rad
  %%pl.t_window = 17.5; % mm/c = ~60 deg of 2.856 GHz
  pl.t_window = 16.7; % mm/c = ~40 deg of 2 GHz (time window to be decided)
  
  window_e1 = pl.e_target*1e3*(1-pl.e_accept);
  window_e2 = pl.e_target*1e3*(1+pl.e_accept);

  % optimise t_ref, e_ref, e_def

  printf("INFO_PL:: optimising t_ref for positron linac . . . \n");

  t_min = floor(min(A_RF(:,5)));
  t_max = ceil(max(A_RF(:,5)));

  np_acc_opt = 0;

  for t_ref = t_min : 1 : t_max
    window_t1 = t_ref - pl.t_window*0.5;
    window_t2 = t_ref + pl.t_window*0.5;
    M_cutT = A_RF(:,5) >= window_t1 & A_RF(:,5) <= window_t2;
    A_RF_cutT = A_RF(M_cutT,:);
    if (rows(A_RF_cutT(:,6))==0) continue; endif
    e_ref = mean(A_RF_cutT(:,6)); 

    % calculate energy deficit

    e_def = 0;
    for i_fix = 1:50
      A_PL_cutT_def = A_RF_cutT;
      A_PL_cutT_def(:,6) = A_RF_cutT(:,6)  + (pl.e_target*1e3 - e_ref + e_def) ...
                           * cos(2 * pi * pl.freq  * (A_RF_cutT(:,5) - t_ref) * 1e-3 / RF_Track.clight);
      M_E_def = A_PL_cutT_def(:,6) >= window_e1 & A_PL_cutT_def(:,6) <= window_e2;
      if (rows(A_PL_cutT_def(M_E_def,6))==0) e_def = e_def + pl.e_target*1e3;
      else e_def = e_def + pl.e_target*1e3 - mean(A_PL_cutT_def(M_E_def,6));
      endif
    endfor

    % limit on energy deficit: 50 MeV

    if (abs(e_def)>50)
      e_def = 50 * e_def/abs(e_def);
    endif

    %% energy calculation with deficit fixed

    A_PL_cutT = A_RF_cutT;
    A_PL_cutT(:,6) = A_RF_cutT(:,6)  + (pl.e_target*1e3 - e_ref + e_def) ...
                     * cos(2 * pi * pl.freq  * (A_RF_cutT(:,5) - t_ref) * 1e-3 / RF_Track.clight);
    M_E = A_PL_cutT(:,6) >= window_e1 & A_PL_cutT(:,6) <= window_e2;
    A_PL_cutTE = A_PL_cutT(M_E,:);
    np_acc = rows(A_PL_cutTE);
    if (np_acc >  np_acc_opt)
      np_acc_opt = np_acc;
      pl.t_ref = t_ref;
      pl.e_ref = e_ref;
      pl.e_def = e_def;
    endif
  endfor

  printf("INFO_PL:: optimised reference particle for positron linac:\n");
  printf("          reference time  : %.1f mm/c\n", pl.t_ref);
  printf("          reference energy: %.1f MeV\n",  pl.e_ref);
  printf("INFO_PL:: energy deficit for positron linac: %.1f MeV\n", pl.e_def);

  %% energy calculation for IL
  window_t1 = pl.t_ref - pl.t_window*0.5;
  window_t2 = pl.t_ref + pl.t_window*0.5;
  A_PL = A_RF;
  A_PL(:,1:4) = 0; % No trans. tracking
  A_PL(:,6) = A_RF(:,6)  + (pl.e_target*1e3 - pl.e_ref + pl.e_def) ...
              * cos(2 * pi * pl.freq  * (A_RF(:,5) - pl.t_ref) * 1e-3 / RF_Track.clight);
  M_cutT = A_PL(:,5) >= window_t1 & A_PL(:,5) <= window_t2;
  M_cutE = A_PL(:,6) >= window_e1 & A_PL(:,6) <= window_e2;
  M_cutTE = M_cutT & M_cutE;
  A_PL_cutT = A_PL(M_cutT,:);
  A_PL_cutE = A_PL(M_cutE,:);
  A_PL_acc = A_PL(M_cutTE,:);
  A_RF_acc = A_RF(M_cutTE,:); % at 200 MeV
  
  % show np
  
  np_pl = rows(A_PL);
  np_pl_cutT = rows(A_PL_cutT);
  np_pl_cutE = rows(A_PL_cutE);
  np_pl_acc  = rows(A_PL_acc);

  printf("INFO:: After IL:\n");
  printf("       MeanE: %.3f [GeV]\n", mean(A_PL(:,6))/1e3);
  printf("       MeanE_cutT: %.3f [GeV]\n", mean(A_PL_cutT(:,6))/1e3);
  printf("       MeanE_cutE: %.3f [GeV]\n", mean(A_PL_cutE(:,6))/1e3);
  printf("       MeanE_acc: %.3f [GeV]\n", mean(A_PL_acc(:,6))/1e3);
  printf("       Accepted E spread: %.1f%%\n", std(A_PL_acc(:,6))/mean(A_PL_acc(:,6))*100);
  printf("       np_pl: %i\n", np_pl);
  printf("       np_pl_cutT: %i\n", np_pl_cutT);
  printf("       np_pl_cutE: %i\n", np_pl_cutE);
  printf("       np_pl_acc: %i\n",np_pl_acc);

  % show accepted emittances
  B6d = Bunch6d(RF_Track.electronmass, 0, +1, A_RF_acc);
  printf("INFO:: Accepted but after RF (200 MeV):\n");
  printf("       Normalised Emitt_X : %.1f mm*rad\n", B6d.get_info().emitt_x/1e3);
  printf("       Normalised Emitt_Y : %.1f mm*rad\n", B6d.get_info().emitt_y/1e3);
  printf("       Normalised Emitt_Z : %.1f mm*rad\n", B6d.get_info().emitt_z/1e3);
  printf("       Normalised Emitt_4D: %.1f mm*rad\n", B6d.get_info().emitt_4d/1e3);
  printf("       Normalised Emitt_6D: %.1f mm*rad\n", B6d.get_info().emitt_6d/1e3);
  printf("       Alpha_X: %.2f\n", B6d.get_info().alpha_x);
  printf("       Alpha_Y: %.2f\n", B6d.get_info().alpha_y);
  printf("       Alpha_Z: %.2f\n", B6d.get_info().alpha_z);
  printf("       Beta_X: %.2f mm/mrad\n", B6d.get_info().beta_x);
  printf("       Beta_Y: %.2f mm/mrad\n", B6d.get_info().beta_y);
  printf("       Beta_Z: %.2f mm/mrad\n", B6d.get_info().beta_z);

  % save

  pl.window_e1 = window_e1;
  pl.window_e2 = window_e2;
  pl.window_t1 = window_t1;
  pl.window_t2 = window_t2;
  pl.np_pl 	= np_pl;
  pl.np_pl_cutT = np_pl_cutT;
  pl.np_pl_cutE = np_pl_cutE;
  pl.np_pl_acc 	= np_pl_acc;

  save('-text',outfname,'A_PL','A_PL_cutT','A_PL_cutE','A_PL_acc','M_cutTE','A_RF_acc','pl');
