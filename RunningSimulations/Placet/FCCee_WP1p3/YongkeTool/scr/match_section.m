
function merit = track_function ( P )

  disp("Free parameters:");
  disp(P);

  global g;

  PM = g.P0;
  PhA = g.PhA;
  Phase = g.Phase;
  Deficit = g.Deficit;

  if (g.isec>=1 && g.isec<=5)
    PM{g.isec} = P;
  elseif (g.isec==6)
    PM = {};
    index = 1;
    for i = 1:columns(g.NPM)
      PM{i} = P(index:index-1+g.NPM(i));
      index = index + g.NPM(i);
    endfor
  elseif (g.isec==7)
    Phase = P(1);
    Deficit = P(2);
  endif

  for i = 1:columns(g.NPM)
    outfile = fopen(["Beamlines/M" int2str(i) ".tcl"], 'w');
    for j = 1:g.NPM(i)
      if (mod(j,2))
        fprintf(outfile, "Drift -length %f\n", PM{i}(j));
      else
        fprintf(outfile, "Quadrupole -length 0.4 -strength %f\n", PM{i}(j));
      endif
    endfor
    fclose(outfile);
  endfor

  outfile = fopen(["Beamlines/free_parameters.tcl"], 'w');
  fprintf(outfile, "set PhA %f\n", PhA);
  fprintf(outfile, "set PHID %f\n", Phase);
  fprintf(outfile, "set DEF_GEV %f\n", Deficit*1e-3);
  fclose(outfile);

  printf("INFO:: Phase: %.1f degree\n", Phase);
  printf("INFO:: Deficit: %.1f MeV\n", Deficit);
  
  %% Distance between quadrupoles > 15 cm and < 5 m
  P_MD = [ PM{1}(1:2:end) PM{2}(1:2:end) PM{3}(1:2:end) PM{4}(1:2:end) PM{5}(1:2:end) ];
  if ( sum(P_MD < 0.15) || sum(P_MD > 5.00) )
    disp("Distance exceeds limits.");
    merit = 1e9 + sum(PM{1}) + sum(PM{2}) + sum(PM{3}) + sum(PM{4}) + sum(PM{5})
    return;
  endif
  
  %% Signs and limits of quadrupole strength
  P_MQF = [ PM{1}(2:4:end) PM{2}(4:4:end) PM{3}(4:4:end) PM{4}(4:4:end) PM{5}(4:4:end) ];
  P_MQD = [ PM{1}(4:4:end) PM{2}(2:4:end) PM{3}(2:4:end) PM{4}(2:4:end) PM{5}(2:4:end) ];
  if ( sum(P_MQF <= 0.01) || sum(P_MQD >= -0.01) )
    disp("Strength exceeds limits.");
    merit = 1e9 + sum(PM{1}) + sum(PM{2}) + sum(PM{3}) + sum(PM{4}) + sum(PM{5})
    return;
  endif

  %% Full tracking
  system(["placet scr/track_section.tcl input/il_input_PLC_Sample.txt " num2str(g.isec)]);

  %% Emittance not Nan value
  E = load("output/Emittance.dat").E;
  if ( sum(isnan(E(:,2)))>0 || sum(isnan(E(:,6)))>0 )
    disp("Emittance has Nan value.");
    merit = 1e9 + sum(PM{1}) + sum(PM{2}) + sum(PM{3}) + sum(PM{4}) + sum(PM{5})
    return;
  endif

  %% Number of particles
  B = load("output/Beam.dat").B;
  NP = rows(B);

  %% Accepted NP
  if (g.isec==7)
    system("octave-cli scr/calc_effective.m");
    NP = load("output/Effective.dat").NP_acc;
  endif

  eff = NP / g.nps;

  merit = -eff

endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option = argv(){1};
IS = argv(){2};

global g;

%% parameters
par_file = [ "par/" option ".m" ];
if(exist(par_file) == 2)
  source(par_file);
endif

g.isec = str2num(IS);

g.NPM = [ columns(g.P_M1) columns(g.P_M2) columns(g.P_M3) columns(g.P_M4) columns(g.P_M5) ];

%% Load MP for matched sections
if (1)
  for i_sec = 1 : g.isec-1
    mp_file = ["output/MP" num2str(i_sec) ".dat"];
    if(exist(mp_file)!=2) continue; endif
    MPi = load(mp_file).MP;
    if(i_sec==1) g.P_M1 = MPi; endif
    if(i_sec==2) g.P_M2 = MPi; endif
    if(i_sec==3) g.P_M3 = MPi; endif
    if(i_sec==4) g.P_M4 = MPi; endif
  endfor
endif

addpath("scr/");

printf("INFO:: Starting point:\n");
printf("       Phase: %.0f degree\n", g.Phase);
printf("       Phase advance: %.0f degree\n", g.PhA);
printf("       Deficit: %.2f\n MeV", g.Deficit);
printf("       Number of matching quadrupoles: %s\n", sprintf("%i ",g.NPM));
printf("       Matching Sec. 1 parameters: %s\n", sprintf("%.4f ",g.P_M1));
printf("       Matching Sec. 2 parameters: %s\n", sprintf("%.4f ",g.P_M2));
printf("       Matching Sec. 3 parameters: %s\n", sprintf("%.4f ",g.P_M3));
printf("       Matching Sec. 4 parameters: %s\n", sprintf("%.4f ",g.P_M4));
printf("       Matching Sec. 5 parameters: %s\n", sprintf("%.4f ",g.P_M5));

g.P0 = { g.P_M1 g.P_M2 g.P_M3 g.P_M4 g.P_M5 };

if(g.isec>=1 && g.isec<=5)
  P = g.P0{g.isec};
elseif(g.isec==6)
  P = [ g.P0{1:end} ];
elseif(g.isec==7)
  P = [ g.Phase g.Deficit ];
else
  P = [];
endif

%% Sample scaling factor
%AI = load("input/il_input_PLC.txt");
AIS = load("input/il_input_PLC_Sample.txt");
g.nps = rows(AIS);
%g.sf = rows(AI) / rows(AIS);
printf("INFO:: Sample NP: %i\n", g.nps);
%printf("INFO:: Sample scaling factor: %.2f\n", g.sf);

options = optimset("Display","iter", "TolX",1e-3, "TolFun",1e-3);

n_loop = 1;
for i_loop = 1:n_loop
  printf("INFO:: Loop %i / %i of fminsearch ...\n", i_loop, n_loop);
  [MP, merit_opt] = fminsearch("track_function", P, options);
  MP = str2num(sprintf('%.4f ', MP));
  disp("Matched parameters:");
  disp(MP);
  P = MP;
endfor

disp("Final matched parameters:");
disp(MP);

save("-text",["output/MP" num2str(g.isec) ".dat"], "MP");
