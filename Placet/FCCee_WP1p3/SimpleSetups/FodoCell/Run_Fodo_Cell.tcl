### Study FODO cells for Linac1 of FCC-ee
###
### Position    betax       alfax   betay       alfay
### start       1.7439 m    0       0.29921 m   0
### end         
### length 1.02 m
### # drift,quad,drift,quad,drift , drift : L in [m], quad : f in [m]
### # last drift fixed by the code to fit total 'length' (see 2 lines above)
### seq 50 -30 30 10


### CONSTANTS
#set pi 3.1415926535897931
#set c 299792458.
#set me 0.0005109989499985808   ;# [GeV]

# FORMULAE FROM THIN LENS APPROXIMATION
# set gammaRel [expr $e_initial / $me]
# set betaRel 0.9999995
# set Rquad 0.1   ;# [m]
# set Bpole 1.   ;# [T]
# set gQuadMax [expr $Bpole / $Rquad]   ;# [T/m]
# set Ra 30.   ;# [mm]
# set Fa 4.
# set muDeg = 76.345   ;# [deg]
# set mu [expr $pi / 180. * muDeg]
# set betax [expr ($Ra/$Fa)**2.*$betaRel*$gammaRel/$emitn_x]   ;# [m]
# set Lcell [expr $beta_x*sin($mu)/(1+sin($mu/2))]   ;# [m]
# set betay [expr $Lcell*(1.-sin($mu/2))/sin($mu)]   ;# [m]
# set fQuad [expr $Lcell/(4.*sin($mu/2))]
# set kQuad [expr $gQuadMax/($e_initial*1e9)*$c]   ;# [1/m2]
# set Lquad [expr 1/($kQuad*$f)]   ;# [m]
# set Ldrift [expr ($Lcell-2*$Lquad)/2]   ;# [m]

### DEFINE INITIAL BEAM PARAMETERS
set charge 1e10   ;# number of particles
set e_initial 0.499489   ;# [GeV]
set e_spread 0.1   ;# [%]
set sigma_z 1.   ;# [um]
set emitn_x 10000.   ;# [mm mrad]
set emitn_y $emitn_x
set beta_x 5.503919   ;# [m]
set beta_y 1.456175   ;# [m]
set alpha_x 0.
set alpha_y 0.
# puts "Initial beam parameters:"
# puts "emitn_x = $emitn_x"
# puts "emitn_y = $emitn_y"
# puts "beta_x = $beta_x"
# puts "beta_y = $beta_y"
# puts "alpha_x = $alpha_x"
# puts "alpha_y = $alpha_y"
# puts "sigma_z = $sigma_z"
# puts "e_spread = $e_spread"
# puts "charge = $charge"

### DEFINE FODO PARAMETERS
set Lquad 1.0   ;# [m]
set fQuad 1.117493   ;# [m]
set Ldrift 0.757640   ;# [m]
# puts "FODO Parameters"
# puts "Lquad = $Lquad m"
# puts "e_initial = $e_initial GeV"
# puts "fQuad = $fQuad m"
# puts "Ldrift = $Ldrift m"

BeamlineNew
SetReferenceEnergy $e_initial
Girder
  Quadrupole -length [expr $Lquad/2.] -strength [expr (+1.)*$e_initial/$fQuad/2.]
  Drift -length $Ldrift
  Quadrupole -length $Lquad -strength [expr (-1.)*$e_initial/$fQuad]
  Drift -length $Ldrift
  Quadrupole -length [expr $Lquad/2.] -strength [expr (+1.)*$e_initial/$fQuad/2.]
BeamlineSet -name beamline


### DEFINE BEAM

array set match {}
set match(emitt_x) [expr $emitn_x*10.]
set match(emitt_y) [expr $emitn_y*10.]
set match(beta_x) $beta_x
set match(beta_y) $beta_y
set match(alpha_x) $alpha_x
set match(alpha_y) $alpha_y
set match(sigma_z) $sigma_z
set match(e_spread) $e_spread
set match(charge) $charge

# Emittance units: [10^-7 m rad]
set common_script_dir ../../scr/common
# set script_dir ../../../../GIT_Placet/examples/fodo_cell
source $common_script_dir/clic_basic_single.tcl
source $common_script_dir/clic_beam.tcl

# Global variables required by function make_beam_many in create_beam.tcl:
# charge, e_initial, match(), n_total
set n_slice 20
set n 1000
set n_total [expr $n_slice * $n]

# Create the beam
make_beam_many beam0 $n_slice $n
FirstOrder 1

# Compute Twiss parameters (method 2)
TwissPlotStep -beam beam0 -file out_twiss_2.dat -step 0.01


### TRACK

Octave {

    # Tracking
    [emitt, beam] = placet_test_no_correction("beamline", "beam0", "None");
    save -text out_beam.dat beam
    save -text out_emitt.dat emitt

    # Compute Twiss parameters (method 1)
    [s, beta_x, beta_y, alpha_x, alpha_y, mu_x, mu_y, Dx, Dy, E] = placet_evolve_beta_function('beamline', $beta_x, $alpha_x, $beta_y, $alpha_y);
    T_1 = [s, beta_x, beta_y, alpha_x, alpha_y, mu_x, mu_y, Dx, Dy, E];
    save -text out_twiss_1.dat T_1

    # Load Twiss parameters from method 2
    T_2 = load('out_twiss_2.dat');

    # Quick analysis
    figure(1)
    subplot(3, 1, 1)
    plot(s, beta_x, 'b-')
    hold on
    plot(s, beta_y, 'r-')
    plot(T_2(:,2), T_2(:,6), 'b--')
    plot(T_2(:,2), T_2(:,10), 'r--')
    xlabel('s [m]')
    ylabel('beta [m]')
    legend('beta_x (Method 1)', 'beta_y (Method 1)', 'beta_x (Method 2)', 'beta_y (Method 2)')
    subplot(3, 1, 2)
    plot(s, alpha_x, 'b-')
    hold on
    plot(s, alpha_y, 'r-')
    plot(T_2(:,2), T_2(:,7), 'b--')
    plot(T_2(:,2), T_2(:,11), 'r--')
    xlabel('s [m]')
    ylabel('alpha')
    legend('alpha_x (Method 1)', 'alpha_y,(Method 1)', 'alpha_x (Method 2)', 'alpha_y (Method 2)')
    subplot(3, 1, 3)
    plot(s, E, 'b-')
    hold on
    plot(T_2(:,2), T_2(:,3), 'b--')
    xlabel('s [m]')
    ylabel('E [GeV]')
    legend('E (Method 1)', 'E (Method 2)')
    #figure(2)
    #plot(beam(:,1), beam(:,2), '.')
    waitforbuttonpress
}
