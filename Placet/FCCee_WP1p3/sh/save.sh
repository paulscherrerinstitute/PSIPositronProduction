#!/bin/bash

## Options
option=CLIC380GeV_Nov2020
#option=CLIC3TeV_Nov2020

[[ $# -eq 1 ]] && option=${1}

#dirname=trk_${option}_final
dirname=${option}

mkdir -pv Results/${dirname}
#rm -rf Results/${dirname}/*
cp -f output/* Results/${dirname}/
