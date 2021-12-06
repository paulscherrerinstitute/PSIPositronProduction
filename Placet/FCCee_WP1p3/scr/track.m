
function merit = track_function ( P )

  global g;

  PM = P;

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
  fprintf(outfile, "set PhA %f\n", g.PhA);
  fprintf(outfile, "set PHID %f\n", g.Phase);
  fprintf(outfile, "set DEF_GEV %f\n", g.Deficit*1e-3);
  fclose(outfile);

  if(g.sampling)
    system(["placet scr/track_section.tcl input/il_input_PLC_Sample.txt"]);
    AI = load("input/il_input_PLC.txt");
    AIS = load("input/il_input_PLC_Sample.txt");
    AOS = load("output/Beam.dat").B;
    sf = rows(AI) / rows(AIS);
    printf("INFO:: Sample scaling factor: %.2f\n", sf);
    system("octave-cli scr/calc_effective.m");
    NP = load("output/Effective.dat").NP_acc;
    EFF_TRANS = NP / rows(AOS)
    EFF_TOT = NP_ACC / rows(AOS)
    NP_ACC = NP * sf
  else
    system("placet scr/track_section.tcl");
  endif

endfunction

%% Tracking

option = argv(){1};

global g;

g.sampling = 0;

%% parameters
par_file = [ "par/" option ".m" ];
if(exist(par_file) == 2)
  source(par_file);
endif

g.NPM = [ columns(g.P_M1) columns(g.P_M2) columns(g.P_M3) columns(g.P_M4) columns(g.P_M5) ];

addpath("scr/");

printf("INFO:: Phase: %.1f degree\n", g.Phase);
printf("INFO:: Phase advance: %.1f degree\n", g.PhA);
printf("INFO:: Deficit: %.2f\n MeV", g.Deficit);
printf("INFO:: Number of matching quadrupoles: %s\n", sprintf("%i ",g.NPM));
printf("INFO:: Matching Sec. 1 parameters: %s\n", sprintf("%.4f ",g.P_M1));
printf("INFO:: Matching Sec. 2 parameters: %s\n", sprintf("%.4f ",g.P_M2));
printf("INFO:: Matching Sec. 3 parameters: %s\n", sprintf("%.4f ",g.P_M3));
printf("INFO:: Matching Sec. 4 parameters: %s\n", sprintf("%.4f ",g.P_M4));
printf("INFO:: Matching Sec. 5 parameters: %s\n", sprintf("%.4f ",g.P_M5));

P0 = { g.P_M1 g.P_M2 g.P_M3 g.P_M4 g.P_M5 };

track_function(P0);
