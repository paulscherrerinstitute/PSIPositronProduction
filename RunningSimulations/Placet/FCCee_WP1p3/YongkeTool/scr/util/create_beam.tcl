## Initialisation

Octave {

  Beta_x = 1;
  Beta_y = 1;
  Alpha_x = 1;
  Alpha_y = 1;
  Emitt_x = 1;
  Emitt_y = 1;
  E_initial = 1;
  E_spread = 1;
  Mean_z = 1;
  Sigma_z = 1;
  Charge = 1;

  printf('create_beam.tcl > Octave_1 > Before RF-Track.\n');
  RF_Track;
  printf('create_beam.tcl > Octave_1 > After RF-Track.\n');

  printf('Loading input_file = $input_file.\n');
  try
    A_PLC = load("$input_file").B;
  catch
    A_PLC = load("$input_file");
  end_try_catch
  A_RFT(:,6) = A_PLC(:,1) * 1e3;  # MeV
  A_RFT(:,1) = A_PLC(:,2) * 1e-3; # mm
  A_RFT(:,3) = A_PLC(:,3) * 1e-3; # mm
  A_RFT(:,5) = A_PLC(:,4) * 1e-3; # mm/c
  A_RFT(:,2) = A_PLC(:,5) * 1e-3; # mrad
  A_RFT(:,4) = A_PLC(:,6) * 1e-3; # mrad

  M = A_RFT(:,5) > -20.000 & A_RFT(:,5) < 50.000;
  A_RFT = A_RFT(M,:);
  disp("MASK")
  size(M)
  sum(M)

  Bunch.mass = RF_Track.electronmass; # MeV/c/c
  Bunch.population = 1e4; # number of particles per bunch
  Bunch.charge     = +1; # units of e+

  % Option 1: Beam input from file
  B_6d = Bunch6d(Bunch.mass, Bunch.population, Bunch.charge, A_RFT);

  % Option 2: Beam input from Twiss parameters
  T = Bunch6d_twiss();
  T.emitt_x = 10000.; % mm.mrad, normalised horizontal emittance
  T.emitt_y = 10000.; % mm.mrad, normalised vertical emittance
  T.emitt_z = 1000.; % mm.permille, normalised longitudinal emittance
  T.alpha_x = 0;
  T.alpha_y = 0;
  T.alpha_z = 0;
  T.beta_x = 1.7439; % m, horizontal beta function
  T.beta_y = 0.29921; % m, vertical beta function
  T.beta_z = 0.1; % m, longitudinal beta function
  T.disp_x = 0; % m, horizontal dispersion
  T.disp_xp = 0; % rad, horizontal dispersion prime
  T.disp_y = 0; % m, vertical dispersion
  T.disp_yp = 0; % rad, vertical dispersion prime
  T.disp_z = 0; % m, longitudinal dispersion
  Pref = 200.; % [MeV/c]
  # Bunch.population is charge in #electrons of the bunch
  Nparticles = 10000;
  B_6d = Bunch6d(Bunch.mass, Bunch.population, Bunch.charge, Pref, T, Nparticles, sigmaCut=0);

  B_6d_info = B_6d.get_info();
  Beta_x = B_6d_info.beta_x; # m
  Beta_y = B_6d_info.beta_y;
  Alpha_x = B_6d_info.alpha_x; #
  Alpha_y = B_6d_info.alpha_y;
  Emitt_x = B_6d_info.emitt_x*1e1; # 1e-7 m.rad
  Emitt_y = B_6d_info.emitt_y*1e1;
  E_initial = B_6d_info.mean_E*1e-3; # GeV
  E_spread = B_6d_info.sigma_E / B_6d_info.mean_E *1e2; # percent
  Mean_z = B_6d_info.mean_t*1e3; # um
  Sigma_z = B_6d_info.sigma_z*1e3; # um
  #Sigma_z = 1000.0*Sigma_z; # um, avoid warning messages (z > 3sigma)
  Charge = size(A_PLC)(1);

  Tcl_SetVar("Beta_x", Beta_x);
  Tcl_SetVar("Beta_y", Beta_y);
  Tcl_SetVar("Alpha_x", Alpha_x);
  Tcl_SetVar("Alpha_y", Alpha_y);
  Tcl_SetVar("Emitt_x", Emitt_x);
  Tcl_SetVar("Emitt_y", Emitt_y);
  Tcl_SetVar("E_initial", E_initial);
  Tcl_SetVar("E_spread", E_spread);
  Tcl_SetVar("Mean_z", Mean_z);
  Tcl_SetVar("Sigma_z", Sigma_z);
  Tcl_SetVar("Charge", Charge);
}

## workaround (tcl bug: can not set inside array)
array set match {
  beta_x 1
  beta_y 1
  alpha_x 1
  alpha_y 1
  emitt_x 1
  emitt_y 1
}
set match(beta_x)  $Beta_x
set match(beta_y)  $Beta_y
set match(alpha_x) $Alpha_x
set match(alpha_y) $Alpha_y
set match(emitt_x) $Emitt_x
set match(emitt_y) $Emitt_y

set e_initial $E_initial ;# GeV
set match(e_spread) $E_spread ;# percent
set match(sigma_z) $Sigma_z ;# um
set match(charge) $Charge ;# number of particles
set charge $match(charge)

puts "INFO:: Mean_z = $Mean_z um."
puts "INFO:: Sigma_z = $Sigma_z um."

source $common_script_dir/clic_basic_single.tcl
source $common_script_dir/clic_beam.tcl

Octave {
  NP = 0;
  NS = 1;
  #infile = fopen("$input_file",'r');
  #NP = fskipl(infile,Inf);
  try
    A = load("$input_file").B;
  catch
    A = load("$input_file");
  end_try_catch
  NP = size(A,1);
  F = factor(NP);
  for i = 1:size(F)
    ns = prod(F(1:i));
    if (prod(F(1:i)) < 20)
      NS = ns;
    else
      break;
    endif
  endfor
  Tcl_SetVar("NP", NP);
  Tcl_SetVar("NS", NS);
}

set n_slice $NS ; # less than 20
set n [ expr $NP / $NS ]
set n_total [expr $n_slice * $n] ; # number of positrons

## Create the beam

make_beam_many beam0 $n_slice $n

FirstOrder 1

## Reset beam (beam0) with input file
Octave {
  printf('create_beam.tcl > Octave_3\n');
  printf('Loading $input_file.\n');
  try
    B0 = load('$input_file').B;
  catch
    B0 = load('$input_file');
  end_try_catch
  disp(B0(1:10,:))
  placet_set_beam("beam0", B0);
}
