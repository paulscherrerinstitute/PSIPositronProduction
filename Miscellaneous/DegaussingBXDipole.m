clear figs

% Parameters for the degaussing cycle
Imax = 200.;
Nstep = 10;
% Computation of the degaussing cycle
n = 0:1:Nstep;
d = 1./(2.^n);
Ideg = Imax * d;
IdegSign = Ideg;
IdegSign(2:2:end) = -Ideg(2:2:end);
% Computation 2 of the degaussing cycle
na=1:1:(Nstep+1);
Idega=Imax*(-1).^(na-1)./(2.^(na-1));

% Expected values from design, assuming linearity
Iexpected = 0:10:200; % [A]
Bexpected = 1.397e3 / 418. * Iexpected; % [mT]
% Measurements during degaussing cycle
ImeasUp1 = [0.5 20:10:100 120:20:200]; % [A]
ImeasDown1 = [180:-20:100 90:-10:0]; % [A]
ImeasUp2 = -1 * [10:10:100]; % [A]
ImeasDown2 = -1 * [90:-10:0]; % [A]
ImeasUp3 = 10:10:50; % [A]
ImeasDown3 = [40:-10:10 5 1 0] % [A]
ImeasUp4 = -1 * [1 5:5:25]; % [A]
ImeasDown4 = -1 * [20:-5:5 1 0]; % [A]
ImeasUp5 = [1 5 10 12.5]; % [A]
ImeasDown5 = [10 5 1 0]; % [A]
ImeasUp6 = -1 * [1 5 6.2]; % [A]
ImeasDown6 = -1 * [5 1 0]; % [A]
ImeasAll = [ ...
  ImeasUp1 ImeasDown1 ImeasUp2 ImeasDown2 ImeasUp3 ImeasDown3 ...
  ImeasUp4 ImeasDown4 ImeasUp5 ImeasDown5 ImeasUp6 ImeasDown6];
BmeasAll = [ ...
  BmeasUp1 BmeasDown1 BmeasUp2 BmeasDown2 BmeasUp3 BmeasDown3 ...
  BmeasUp4 BmeasDown4 BmeasUp5 BmeasDown5 BmeasUp6 BmeasDown6];
BmeasUp1 = [0.25 6.1 9 12 15 17 21 23 27 30 36 42.9 48.0 54.0 59.9]; % [mT]
BmeasDown1 = [54.1 48.1 42.1 36.2 30.2 27.2 24.2 21.2 18.1 15.2 12.2 9.24 6.19 3.18 0.12]; % [mT]
BmeasUp2 = -1 * [2.97 6.07 9.12 12.17 15.20 18.3 21.4 24.5 27.5 30.6]; % [mT]
BmeasDown2 = -1 * [27.6 24.5 21.5 18.4 15.4 12.31 9.24 6.20 3.20 0.12]; % [mT]
BmeasUp3 = [3.02 6.04 9.07 12.11 15.12]; % [mT]
BmeasDown3 = [12.13 9.14 6.12 3.10 1.59 0.38 0.04]; % [mT]
BmeasUp4 = -1 * [0.31 1.52 3.05 4.56 6.08 7.61]; % [mT]
BmeasDown4 = -1 * [6.12 4.61 3.11 1.61 0.38 0.07]; % [mT]
BmeasUp5 = [0.28 1.50 3.01 3.75]; % [mT]
BmeasDown5 = [3.01 1.54 0.32 0.00]; % [mT]
BmeasUp6 = -1 * [0.34 1.54 1.90]; % [mT]
BmeasDown6 = -1 * [1.56 0.36 0.02]; % [mT]

% Plot degaussing cycle
##figure(1)
##plot(n, d, 'o')
figure(2)
plot(n+1, Ideg, 'o')
hold on
plot(n+1, IdegSign, '*')
plot(na, Idega, '.g')
hold off
legend('Ideg', 'IdegSign', 'Idega')
xlabel('Step number')
ylabel('Current [A]')
grid
##set(gca, 'FontSize', 24)

fprintf('%.1f ', IdegSign')
fprintf('\n')

% Plot performed ramps
figure(3)
plot(ImeasUp1, BmeasUp1, '.-')
hold on
plot(ImeasDown1, BmeasDown1, '.-')
plot(-ImeasUp2, -BmeasUp2, '.--')
plot(-ImeasDown2, -BmeasDown2, '.--')
plot(ImeasUp3, BmeasUp3, '.-')
plot(ImeasDown3, BmeasDown3, '.-')
plot(-ImeasUp4, -BmeasUp4, '.--')
plot(-ImeasDown4, -BmeasDown4, '.--')
plot(ImeasUp5, BmeasUp5, '.-')
plot(ImeasDown5, BmeasDown5, '.-')
plot(-ImeasUp6, -BmeasUp6, '.--')
plot(-ImeasDown6, -BmeasDown6, '.--')
plot(Iexpected, Bexpected/11, 'o:')
hold off
xlabel('Current [A]')
ylabel('Magnetic field [mT]')
grid
##set(gca, 'FontSize', 24)

fprintf('%20s %20s\n', 'Imeas [A]', 'Bmeas [mT]')
fprintf('%20.1f %20.1f\n', [ImeasAll; BmeasAll])
fprintf('\n')

