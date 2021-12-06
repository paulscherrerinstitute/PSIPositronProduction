#!/bin/bash

option=CLIC380GeV_Nov2020
[[ $# -eq 1 ]] && option=${1}

###################################################
## {1..5}: Sec 1-5;  6: All Secs;  7: Phase Deficit
###################################################
for IS in 7
do
  [[ $IS -ge 1 && $IS -le 5 ]] && echo Doing matching for Section $IS ..
  [[ $IS -eq 6 ]] && echo Doing matching for all Sections ..
  [[ $IS -eq 7 ]] && echo Optimising Phase and Deficit ...
  octave-cli scr/match_section.m $option $IS |& tee output/log.match_section
done
