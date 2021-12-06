
%option = "CLIC380GeV_Nov2020";
%option = "CLIC3TeV_Nov2020";
%option = "CLIC380GeV_Nov2020_HugoLinAMD";
%option = "CLIC3TeV_Nov2020_HugoLinAMD";
%option = "CLIC380GeV_Nov2020_HugoAMD";
option = "CLIC3TeV_Nov2020_HugoAMD";

filename = ["../IteScan/job/Dat/TW_trk_" option "_final_0_0.dat"];

A = load(filename).A_TW;

disp("Mean(A)")
mean(A)
disp("Std(A)")
std(A)

if(1)
  %hist(A(:,6), [0+25:50:500+25]);
  %hist(A(:,6));
  %hist(A(:,6), [0:50:500]); %% E
  %hist(A(:,2), [-14:2:14]); %% xp
  hist(A(:,4), [-16:2:16]); %% xp
endif

if(0)
  %scatter(A(:,5), A(:,6));
endif

if(0)
  %scatter(A(:,5), A(:,6));
endif

pause();
