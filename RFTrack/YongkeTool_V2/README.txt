
1) positron distributions:

  distributions/positrons_200MeV_CLIC_Lband.dat is for HTS + CLIC L-band @ 0.5 T
  distributions/positrons_200MeV_PSI_Sband.dat is for HTS + PSI S-band @ 1.5 T

  Output is in Octave / Matlab text format.
  
  A_RF and A_RF_ACC are bunches at the end of capture RF linac. 
  A_RF contains all positrons. 
  A_RF_ACC contains positrons that are finally accepted by DR.
  
  Positrons are saved in different rows.
  Columns are: x/mm x'/mrad y/mm y'/mrad t/(mm/c) p/(MeV/c)


2) Codes:

  track_CLIC_Lband_200MeV.m is for tracking of CLIC L-band
  track_PSI_Sband_200MeV.m is for tracking of PSI S-band

3) Field maps:

  field/field_map_HTS_5coils.dat is HTS's field map
  field/field_map_CLIC_Lband.dat is CLIC L-band's field map
  field/field_map_PSI_Sband.dat is PSI S-band's field map
