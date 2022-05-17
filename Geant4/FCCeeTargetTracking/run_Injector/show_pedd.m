amor_file = "/afs/psi.ch/project/Pcubed/SimulationRuns/Geant4/000007/FCCeeTargetTracking_amor.dat";

A_amor = load(amor_file);

Ne_s = 10000; %% number of simulated e-

X0 = 3.5; %% mm, radiation length
density = 19.25;  %% g/cm^3, W density
amor_thick = 17.5; %% mm

nbins_xy = round(amor_thick*2 + 1);
delta_xy = (amor_thick + 0.5) / nbins_xy;
delta_x = delta_xy;
delta_y = delta_xy;

nbins_z = round(amor_thick*2 + 1);
delta_z = (amor_thick + 0.5) / nbins_z;
volume_cell_amor = delta_x * delta_y * delta_z * 1.e-3; 	%% cm^3

factor_GeV_to_J = 1.60218e-10;
factor_GeV_to_kW = factor_GeV_to_J*1.0e-3; %% per second

[peak_energy_amor, i_max_amor] = max(A_amor(:,4) * 1e-3); %% GeV

total_energy_amor = sum(A_amor(:,4) * 1e-3); %% GeV

yield = 2.0; %% assume accepted yield is 1.
Qe = 1.602e-19;   % Electron charge (without sign)

% Case FCC-ee
%Nb_e  = 2;   % Number of e- bunches per pulse
%Np_b  = 3.4e-9 / Qe;    % Required number of e+ per bunch
%Ne_b  = Np_b / yield;   % Number of e- per bunch
%Ne_b  = 1.19e-9 / Qe;    % Number of e- per bunch
%f_rep = 200;   % Repetition rate [Hz]

% Case P3
Nb_e  = 1;   % Number of e- bunches per pulse
Ne_b  = 200e-12 / Qe;   % Number of e- per bunch
f_rep = 1;   % Repetition rate [Hz]

pedd_amor = peak_energy_amor * (Nb_e*Ne_b/Ne_s) / (volume_cell_amor*density) * factor_GeV_to_J;
power_amor = total_energy_amor * (f_rep*Nb_e*Ne_b/Ne_s) * factor_GeV_to_kW;

printf("Assuming final accepted yield by DR is %.2f:\n",yield);
printf("  Amorphous target PEDD: %f [J/g]\n", pedd_amor);
printf("  Amorphous target total deposited power: %f [kW]\n", power_amor);

