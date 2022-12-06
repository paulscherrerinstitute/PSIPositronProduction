
%A = load("output/Emittance.dat").E;
A = load("Results/CLIC380GeV_Nov2020/Emittance.dat").E;
%A = load("Results/CLIC380GeV_Nov2020_HugoLinAMD/Emittance.dat").E;

for i = 1 : size(A)(1)
  A(i,2) = A(i,2) * 1e-4; %% 10^-7 m.rad -> mm.rad
  A(i,6) = A(i,6) * 1e-4;
endfor

exmax = 1e+6;
eymax = 1e+6;

if(A(1,1)==0)
  A(1,1) = 1e-3;
endif

if (max(A(:,2)) < exmax)
  exmax = max(A(:,2));
endif

if (max(A(:,6)) < eymax)
  eymax = max(A(:,6));
endif

fig = figure();
%fig = figure("visible","off"); % display off

subplot(2,2,1);
%scatter(A(:,1),A(:,2));
plot(A(:,1),A(:,2));
set(gca, "xscale","lin");
ylim([0, 1.2*exmax]);
ylabel("X emittance [mm.rad]");
xlabel("z [m]");

subplot(2,2,2);
%scatter(A(:,1),A(:,6));
plot(A(:,1),A(:,6));
set(gca, "xscale","lin");
ylim([0, 1.2*eymax]);
ylabel("Y emittance [mm.rad]");
xlabel("z [m]");

subplot(2,2,3);
%scatter(A(:,1),A(:,2));
plot(A(:,1),A(:,2));
set(gca, "xscale","log");
ylim([0, 1.2*exmax]);
ylabel("X emittance [mm.rad]");
xlabel("z [m]");

subplot(2,2,4);
%scatter(A(:,1),A(:,6));
plot(A(:,1),A(:,6));
set(gca, "xscale","log");
ylim([0, 1.2*eymax]);
ylabel("Y emittance [mm.rad]");
xlabel("z [m]");

%% Pause

pause();

%% Save

%print(fig,"output/figure.pdf");
%print(fig,"output/emittance_evolution.png");
%open output/figure.pdf
