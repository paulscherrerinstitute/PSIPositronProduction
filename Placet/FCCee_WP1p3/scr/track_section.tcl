## Number of CPU cores to use

ParallelThreads -num 32

## input beam
if {$argc > 0} {
  set input_file [lindex $argv 0]
} else {
  set input_file "input/il_input_PLC.txt"
}

## section number
if {$argc > 1} {
  set i_sec [lindex $argv 1]
} else {
  set i_sec 0
}

## Path variables

set common_script_dir scr/common
set util_script_dir scr/util
set beamline_dir Beamlines

############ Start of beamline ##################

source $beamline_dir/free_parameters.tcl
source $beamline_dir/parameters.tcl

BeamlineNew
Girder
Quadrupole -length 0 -strength 0

if {$i_sec ==0 || $i_sec >= 1} {
  source $beamline_dir/M1.tcl
  source $util_script_dir/create_S1.tcl
}
if {$i_sec ==0 || $i_sec >= 2} {
  source $beamline_dir/M2.tcl
  source $util_script_dir/create_S2.tcl
}
if {$i_sec ==0 || $i_sec >= 3} {
  source $beamline_dir/M3.tcl
  source $util_script_dir/create_S3.tcl
}
if {$i_sec ==0 || $i_sec >= 4} {
  source $beamline_dir/M4.tcl
  source $util_script_dir/create_S4.tcl
}
if {$i_sec ==0 || $i_sec >= 5} {
  source $beamline_dir/M5.tcl
  source $util_script_dir/create_S5.tcl
}

BeamlineSet -name beamline

############ End of beamline ##################

## Set aperture

source $util_script_dir/set_aperture.tcl

## Create beam

source $util_script_dir/create_beam.tcl

## Tracking & Save

Octave {

  output_beam_format = "%s %ex %sex %x %sx %ey %sey %Env %sy %n";

  [E, B] = placet_test_no_correction("beamline", "beam0", "None",output_beam_format);

  save -text output/Emittance.dat E;
  save -text output/Beam.dat B;

  ## tracking only
  if ($i_sec == 0)
    save -text output/Emittance.dat E;
    save -text output/Beam.dat B;

    T = placet_get_twiss_matrix(B);
    save -text output/Twiss_All.dat T;

    I = placet_get_name_number_list('beamline','*');
    S = placet_element_get_attribute('beamline', I, 's'); # m
    LOSS = placet_element_get_attribute('beamline', I, 'aperture_losses'); # perc
    A_LOSS = [ S LOSS ];
    save -text output/Loss.dat A_LOSS;

    ## what's this?
    #TwissPlot -beam beam0 -file output/twiss_plot.dat
  endif
}
