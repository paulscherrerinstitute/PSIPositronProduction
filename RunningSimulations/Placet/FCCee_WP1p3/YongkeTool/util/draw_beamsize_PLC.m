
%A = load("output/Emittance.dat").E;
A = load("Results/CLIC380GeV_Nov2020/Emittance.dat").E;
%A = load("Results/CLIC380GeV_Nov2020_HugoLinAMD/Emittance.dat").E;

jx = 5;
jy = 9;

for i = 1 : size(A)(1)
  A(i,jx) = A(i,jx) * 1e-4; %% um -> cm
  A(i,jy) = A(i,jy) * 1e-4;
endfor

exmax = 1.0; %% cm
eymax = 1.0;

if (max(A(:,jx)) > exmax)
  exmax = max(A(:,jx));
endif

if (max(A(:,jy)) > eymax)
  eymax = max(A(:,jy));
endif

fig = figure();
%fig = figure("visible","off"); % display off

subplot(2,2,1);
%scatter(A(:,1),A(:,jx));
plot(A(:,1),A(:,jx));
set(gca, "xscale","lin");
ylim([0, 1.2*exmax]);
ylabel("X beam size [cm]");
xlabel("z [m]");

subplot(2,2,2);
%scatter(A(:,1),A(:,jy));
plot(A(:,1),A(:,jy));
set(gca, "xscale","lin");
ylim([0, 1.2*eymax]);
ylabel("Y beam size [cm]");
xlabel("z [m]");

subplot(2,2,3);
%scatter(A(:,1),A(:,jx));
plot(A(:,1),A(:,jx));
set(gca, "xscale","log");
ylim([0, 1.2*exmax]);
ylabel("X beam size [cm]");
xlabel("z [m]");

subplot(2,2,4);
%scatter(A(:,1),A(:,jy));
plot(A(:,1),A(:,jy));
set(gca, "xscale","log");
ylim([0, 1.2*eymax]);
ylabel("Y beam size [cm]");
xlabel("z [m]");

%% Pause

pause();

%% Save

%print(fig,"output/beamsize_evolution.png");
%open output/figure.pdf
