

B = load('distributions/positrons_200MeV_CLIC_Lband.dat').A_RF;

%

%X  = B(:,1); % mm
%XP = B(:,2); % mrad
%Y  = B(:,3); % mm
%YP = B(:,4); % mrad
T  = B(:,5); % mm/c
P  = B(:,6); % MeV/c

% cuts
M = P < 350 & T > 18785 & T < 19375;
T  = T(M,:);
P  = P(M,:);

figure('visible','off');

nr = 2;
nc = 1;
ip = 1;

pkg load statistics;

% draw

subplot(nr,nc, ip++);
  scatter(T,P);
  ylabel('P [MeV/c]');
  xlabel('T [mm/c]');
  
subplot(nr,nc, ip++);
  [NN, CC] = hist3 ([T P], [200 200]);
  imagesc (CC{1}, CC{2}, NN');
  axis xy;
  colorbar;
  ylabel('P [MeV/c]');
  xlabel('T [mm/c]');


% plot

print('plot.png', '-dpng', '-S1600,900');

system('feh plot.png');
system('rm plot.png');

