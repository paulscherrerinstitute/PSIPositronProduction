#!/bin/bash

## options
Options=(
  #CLIC380GeV_Nov2020
  CLIC3TeV_Nov2020
  #CLIC380GeV_Nov2020_HugoLinAMD
  #CLIC3TeV_Nov2020_HugoLinAMD
  #CLIC380GeV_Nov2020_HugoAMD
  #CLIC3TeV_Nov2020_HugoAMD
)

do_match=0

for option in ${Options[*]}
do

  rm -rf output/*
  ./sh/prepare.sh $option
  
  if [[ $do_match -eq 1 ]];then
    ./sh/match.sh $option
  else
    ./sh/track.sh $option
  fi
  
  ./sh/save.sh $option

done
