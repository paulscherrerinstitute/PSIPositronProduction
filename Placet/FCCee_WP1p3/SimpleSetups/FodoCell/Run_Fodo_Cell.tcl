### Study FODO cells for Linac1 of FCC-ee
###
### Position    betax       alfax   betay       alfay
### start       1.7439 m    0       0.29921 m   0
### end         
### length 1.02 m
### # drift,quad,drift,quad,drift , drift : L in [m], quad : f in [m]
### # last drift fixed by the code to fit total 'length' (see 2 lines above)
### seq 50 -30 30 10


### DEFINE BEAMLINE

set charge 10000   ;# number of particles
set e_initial 0.2   ;# [GeV]
set gammaRel 392.390
set betaRel 1.
set e_spread 0.1   ; # [%]
set sigma_z 1.   ;# [um]
set emitnx 10000.   ;# [mm mrad]
set emitny $emitnx
set alphax 0.
set alphay 0.
set Ra 20.   ;# [mm]
set Fa 3.
set betax [expr ($Ra/$Fa)**2.*$betaRel*$gammaRel/$emitnx]   ;# [m]
set Lcell [expr $betax/(1+sqrt(2)/2)]   ;# [m]
set betay [expr $Lcell*(1-sqrt(2)/2)]   ;# [m]
set f [expr $Lcell/(4.*1.)]
set gQuad 6.0   ;# [T/m]
set Lquad 0.30776   ;# [m]
set kQuad 8.9964   ;# [1/m]
set Ldrift [expr ($Lcell-2*$Lquad)/2]   ;# [m]

# Parameters of original example
# set charge 2.56e9
# set e_initial 9.0
# set e_spread 2.0
# set sigma_z 0
# set betax 170.71067812
# set betay 29.28932188 
# set alphax 2.41421356  
# set alphay -0.41421356
# set Lquad 1.0
# set kQuad 30.0
# set Ldrift 30.0
# set emitnx 0.1
# set emitny 0.1

#puts "Lcell = $Lcell m"
puts "Lquad = $Lquad m"
puts "Ldrift = $Ldrift m"
puts "betax = $betax m"
puts "betay = $betay m"

BeamlineNew
Girder
  Quadrupole -length [expr $Lquad/2.] -strength [expr (+1)*$e_initial*$kQuad*$Lquad/2.]
  Drift -length $Ldrift
  Quadrupole -length $Lquad -strength [expr (-1.)*$e_initial*$kQuad*$Lquad]
  Drift -length $Ldrift
  Quadrupole -length [expr $Lquad/2.] -strength [expr (+1)*$e_initial*$kQuad*$Lquad/2.]
BeamlineSet -name beamline


### DEFINE BEAM

array set match {}
set match(emitt_x) [expr $emitnx*10.]
set match(emitt_y) [expr $emitny*10.]
set match(beta_x) $betax
set match(beta_y) $betay
set match(alpha_x) $alphax
set match(alpha_y) $alphay
set match(sigma_z) $sigma_z
set match(e_spread) $e_spread
set match(charge) $charge
# puts $match(emitt_x)
# puts $match(emitt_y)
# puts $match(beta_x)
# puts $match(beta_y)
# puts $match(alpha_x)
# puts $match(alpha_y)
# puts $match(sigma_z)
# puts $match(e_spread)
# puts $match(charge)

# Emittance units: [10^-7 m rad]
set common_script_dir ../scr/common
# set script_dir ../../../../GIT_Placet/examples/fodo_cell
source $common_script_dir/clic_basic_single.tcl
source $common_script_dir/clic_beam.tcl

# Global variables required by function make_beam_many in create_beam.tcl:
# charge, e_initial, match(), n_total
set n_slice 1
set n $charge
set n_total [expr $n_slice * $n]

# Create the beam
make_beam_many beam0 $n_slice $n
FirstOrder 1


### TRACK

Octave {

    # Tracking
    [emitt,beam] = placet_test_no_correction("beamline", "beam0", "None");
    [s, beta_x, beta_y, alpha_x, alpha_y, mu_x, mu_y, Dx, Dy, E] = placet_evolve_beta_function('beamline', $betax, $alphax, $betay, $alphay);

    # Save output
    T = [s, beta_x, beta_y, alpha_x, alpha_y, mu_x, mu_y, Dx, Dy, E];
    save -text output_twiss.dat T;
    save -text beam.out beam

    # Quick analysis
    disp(size(beam))
    figure(1)
    subplot(2, 1, 1)
    plot(s, beta_x)
    hold on
    plot(s, beta_y)
    subplot(2, 1, 2)
    plot(s, alpha_x)
    hold on
    plot(s, alpha_y)
    figure(2)
    plot(beam(:,1), beam(:,2), '.')
    waitforbuttonpress
}
