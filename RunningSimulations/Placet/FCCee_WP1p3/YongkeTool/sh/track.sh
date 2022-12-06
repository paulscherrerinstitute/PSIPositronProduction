#!/bin/bash

option=CLIC380GeV_Nov2020
[[ $# -eq 1 ]] && option=${1}

octave-cli scr/track.m $option |& tee output/log.track

octave-cli scr/calc_effective.m |& tee output/log.calc_effective
