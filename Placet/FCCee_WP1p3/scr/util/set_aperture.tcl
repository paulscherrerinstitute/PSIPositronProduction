Octave {
  RQ = 3e-2; # m
  RC = 2e-2; # m
  RD = 3e-2; # m
  NL_Q = placet_get_number_list('beamline', 'quadrupole');
  NL_C = placet_get_number_list('beamline', 'cavity');
  NL_D = placet_get_number_list('beamline', 'drift');
  placet_element_set_attribute('beamline', NL_Q, 'aperture_x', RQ); # m
  placet_element_set_attribute('beamline', NL_Q, 'aperture_y', RQ);
  placet_element_set_attribute('beamline', NL_C, 'aperture_x', RC);
  placet_element_set_attribute('beamline', NL_C, 'aperture_y', RC);
  placet_element_set_attribute('beamline', NL_D, 'aperture_x', RD);
  placet_element_set_attribute('beamline', NL_D, 'aperture_y', RD);
  placet_element_set_attribute('beamline', NL_Q, 'aperture_shape', 'circular');
  placet_element_set_attribute('beamline', NL_C, 'aperture_shape', 'circular');
  placet_element_set_attribute('beamline', NL_D, 'aperture_shape', 'circular');
}
