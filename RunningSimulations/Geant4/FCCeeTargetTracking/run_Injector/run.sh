#!/bin/bash

config_file=config.mac
output_filename=FCCeeTargetTracking

ncore=$(nproc --all)
  
## options: "all", "primary", "photon_emit", "xtal_leave", "amor_arrive", "amor_leave", "amd_arrive","amd_leave"
tree_option="all"
seed=1
  
../Injector_build/injector $config_file ${ncore} $tree_option $seed |& tee logfile
  
## Merge outputs
echo "Merging root files ..."
if [[ $ncore -gt 1 ]];then
  hadd -f ${output_filename}.root ${output_filename}_t*.root
  rm -f ${output_filename}_t*.root
elif [[ $ncore -eq 1 ]];then
  mv ${output_filename}_t0.root ${output_filename}.root
fi
 
root -l -b -q show_N_positrons.C

## Convert to Pcubed standard format
python convert_fcceett_to_standard_df.py ${output_filename}.root
