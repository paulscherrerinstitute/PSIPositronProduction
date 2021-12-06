#!/bin/bash

## options
option=CLIC380GeV_Nov2020
#option=CLIC3TeV_Nov2020

[[ $# -eq 1 ]] && option=${1}

suffix=trk_${option}_final
input_RFT=../IteScan/job/Dat/TW_${suffix}_0_0.dat

octave-cli scr/prepare_input.m $input_RFT |& tee output/log.prepare_input

#cp -f ../IteScan/job/Results/${suffix}_0_0.dat input/analytic_result.dat
