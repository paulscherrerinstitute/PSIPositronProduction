
A = load(["output/Beam.dat"]).B;

NP_out = size(A)(1);

A_in = load("input/il_input_PLC.txt");
NP_in = size(A_in)(1);

%% energy acceptance

A_E = A(:,1);
M_E = abs(A_E/2.86 - 1.0) < 1.2e-2;
printf("INFO:: E window: [%f, %f] GeV.\n",2.86*(1.0-1.2e-2),2.86*(1.0+1.2e-2));

%% aperture

A_X = A(:,2);
A_Y = A(:,3);
M_R = sqrt(A_X.*A_X+A_Y.*A_Y)*1e-3 < 20;

%%% time / z window

%% find time window

if(1)

  M = M_E & M_R;
  A_Z = A(M,4) * 1e-3;

  zmin = min(A_Z);
  zmax = max(A_Z);

  if (size(A_Z)(1) == 0)
    disp("Accepted positrons: 0");
    return;
  endif

  z1 = zmin;
  dZ = 19.8; %% mm

  z1_opt = z1;
  np_opt = 0;

  for z1 = zmin : 0.1 : zmax
  
   M_Z = A_Z > z1 & A_Z < z1 + dZ;
   A_Z_eff = A_Z(M_Z);

   np = size(A_Z_eff)(1);
  
   if (np>np_opt)
     np_opt = np;
     z1_opt = z1;
   endif
  endfor

  printf("INFO:: z1_opt: %f mm.\n", z1_opt);
  printf("INFO:: Z window: [%f, %f] mm.\n", z1_opt, z1_opt+dZ);

  A_Z = A(:,4);
  M_Z = A_Z * 1e-3 > z1_opt & A_Z * 1e-3 < z1_opt + dZ;

endif

M = M_E & M_Z & M_R;

A_eff = A(M,:);


NP_acc = size(A_eff)(1);

printf("INFO:: Input NP = %i\n",NP_in);
printf("INFO:: Output NP = %i\n",NP_out);
printf("INFO:: Accepted NP = %i\n",NP_acc);
printf("INFO:: Transport efficiency = %f\n",(1.0*NP_out/NP_in));
printf("INFO:: Acceptance efficiency = %f\n",(1.0*NP_acc/NP_out));
printf("INFO:: Total efficiency = %f\n",(1.0*NP_acc/NP_in));

min_e_eff = min(A_eff(:,1));
max_e_eff = max(A_eff(:,1));
min_z_eff = min(A_eff(:,4))*1e-3;
max_z_eff = max(A_eff(:,4))*1e-3;
mean_e_eff = mean(A_eff(:,1));
printf("INFO:: min_e (accepted) = %f GeV\n",min_e_eff);
printf("INFO:: max_e (accepted) = %f GeV\n",max_e_eff);
printf("INFO:: min_z (accepted) = %f mm\n",min_e_eff);
printf("INFO:: max_z (accepted) = %f mm\n",max_z_eff);
printf("INFO:: mean_e (accepted) = %f GeV\n",mean_e_eff);

save("output/Effective.dat","z1_opt","NP_in","NP_out","NP_acc","min_e_eff","max_e_eff","min_z_eff","max_z_eff","mean_e_eff");
