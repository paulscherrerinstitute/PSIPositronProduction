
input_file = argv(){1};

printf("INFO:: Input file: %s.\n", input_file);

output_file = "input/il_input_PLC.txt";

A_TW_origin = load(input_file).A_TW;
printf("INFO:: Original input beam:\n");
printf("       Number of e+ is %i \n",size(A_TW_origin)(1));
printf("       Mean energy: %f MeV\n", mean(A_TW_origin(:,6)) );
printf("       Energy spread: %f\n", (std(A_TW_origin(:,6))/mean(A_TW_origin(:,6))) );

A_TW_window = A_TW_origin;

%% Apply 3 sigma window
if(1)
  Mean = mean(A_TW_origin);
  Sigma = std(A_TW_origin);
  M_E = abs(A_TW_origin(:,6) - Mean(6)) < 3.0*Sigma(6);
  M_XP = abs(A_TW_origin(:,2) - Mean(2)) < 3.0*Sigma(2);
  M_YP = abs(A_TW_origin(:,4) - Mean(4)) < 3.0*Sigma(4);
  M = M_E & M_XP & M_YP;
  A_TW_window = A_TW_origin(M,:);
  printf("INFO:: Input beam with 3 sigma [e, xp, yp] window applied:\n");
  printf("       Number of e+ is %i \n",size(A_TW_window)(1));
  printf("       Mean energy: %f MeV\n", mean(A_TW_window(:,6)) );
  printf("       Energy spread: %f\n", (std(A_TW_window(:,6))/mean(A_TW_window(:,6))) );
endif

%% Apply energy window
if(0)
  e_min = 30; %% MeV
  e_max = 400; %% MeV
  M_E = A_TW_origin(:,6) > e_min & A_TW_origin(:,6) < e_max; %% MeV
  M = M_E;
  A_TW_window = A_TW_origin(M,:);
  printf("INFO:: Input beam with [%i, %i] MeV energy window applied:\n", e_min, e_max);
  printf("       Number of e+ is %i \n",size(A_TW_window)(1));
  printf("       Mean energy: %f MeV\n", mean(A_TW_window(:,6)) );
  printf("       Energy spread: %f\n", (std(A_TW_window(:,6))/mean(A_TW_window(:,6))) );
endif

%% Apply time window
if(0)
  %Tarr = 18733; %% mm/c (opt3TeV5GeV2019)
  Tarr = 18708; %% mm/c (opt3TeV5GeV2020)
  M_T = A_TW_window(:,5) >= Tarr & A_TW_window(:,5) <= Tarr + 20;
  M = M_T;
  A_TW_window = A_TW_window(M,:);
  printf("INFO:: Input beam with 20 mm/c time window applied:\n");
  printf("       Number of e+ is %i \n",size(A_TW_window)(1));
  printf("       Mean energy: %f MeV\n", mean(A_TW_window(:,6)) );
  printf("       Energy spread: %f\n", (std(A_TW_window(:,6))/mean(A_TW_window(:,6))) );
endif

%% Shrink energy spread
if(0)
  A_TW_window_shrink = A_TW_window;
  for i = 1:size(A_TW_window)(1)
    A_TW_window_shrink(i,6) = mean(A_TW_window(:,6)) + 0.30*( A_TW_window(i,6) - mean(A_TW_window(:,6)) );
  endfor
  A_TW_window = A_TW_window_shrink;
  printf("INFO:: Input beam with energy spread shrinked:\n");
  printf("       Number of e+ is %i \n",size(A_TW_window)(1));
  printf("       Mean energy: %f MeV\n", mean(A_TW_window(:,6)) );
  printf("       Energy spread: %f\n", (std(A_TW_window(:,6))/mean(A_TW_window(:,6))) );
endif

RF_Track;

Bunch.mass = RF_Track.electronmass; % MeV/c/c
g.population = 1e4;
Bunch.population = g.population; %% number of particles per bunch
Bunch.charge     = +1; %% units of e+
      
B_TW = Bunch6d(Bunch.mass, Bunch.population, Bunch.charge, A_TW_window);

%% Input information
B_TW_info = B_TW.get_info();
Beta_x = B_TW_info.beta_x; % m
Beta_y = B_TW_info.beta_y;
Alpha_x = B_TW_info.alpha_x; %
Alpha_y = B_TW_info.alpha_y;
Emitt_x = B_TW_info.emitt_x*1e1; % 1e-7 m.rad
Emitt_y = B_TW_info.emitt_y*1e1;
E_initial = B_TW_info.mean_E*1e-3; % GeV
E_spread = B_TW_info.sigma_E / B_TW_info.mean_E *1e2; % percent
Mean_z = B_TW_info.mean_t*1e3; % um
Sigma_z = B_TW_info.sigma_z*1e3; % um
printf("INFO:: Input information:\n");
printf("       Beta_x: %f [m]\n", Beta_x);
printf("       Beta_y: %f [m]\n", Beta_y);
printf("       Alpha_x: %f [rad]\n", Alpha_x);
printf("       Alpha_y: %f [rad]\n", Alpha_y);
printf("       Emitt_x: %f [1e-7.m.rad]\n", Emitt_x);
printf("       Emitt_y: %f [1e-7.m.rad]\n", Emitt_y);
printf("       Mean_z: %f [um]\n", Mean_z);
printf("       Sigma_z: %f [um]\n", Sigma_z);

%% Placet format
A_TW = B_TW.get_phase_space("%Pc %x %y %t %xp %yp");
A_TW_placet(:,1) = A_TW(:,1) * 1e-3; %% MeV to GeV
A_TW_placet(:,2) = A_TW(:,2) * 1e3; %% mm to um
A_TW_placet(:,3) = A_TW(:,3) * 1e3; 
A_TW_placet(:,4) = (A_TW(:,4) - mean(A_TW(:,4)))* 1e3; %% mm/c to um
A_TW_placet(:,5) = A_TW(:,5) * 1e3; %% mrad to urad
A_TW_placet(:,6) = A_TW(:,6) * 1e3;

%% Input for tracking
csvwrite(output_file, A_TW_placet)

%% Input for matching
A = load(output_file);
pkg load statistics;
na = rows(A);
ma = 3000;
IS = randsample(na, ma);
AS = A(IS,:);
output_file = "input/il_input_PLC_Sample.txt";
csvwrite(output_file,AS);
printf("INFO:: Sampling for matching:\n");
printf("       Sample size: %i\n", ma);
