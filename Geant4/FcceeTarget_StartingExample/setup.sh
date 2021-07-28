#!/bin/bash

## on lxplus ##

if [[ $(hostname) == "lxplus6"* ]]; then

  export WORK_DIR=/afs/cern.ch/work/y/yozhao

  ## CLIC env ##
  
  ## setup latest GCC, CMake, ROOT, Geant4, Octave, Etc.
  
  source /cvmfs/clicbp.cern.ch/x86_64-slc6-gcc8-opt/setup.sh
  
  ## Variables for Geant4
  
  export G4_CMAKE_DIR=$G4INSTALL/lib64/Geant4-$(geant4-config --version)
  
  ## To make condor work  without errors ##
  
  unset PYTHONPATH PYTHONHOME

elif [[ $(hostname) == "lxplus7"* ]]; then

  export WORK_DIR=/afs/cern.ch/work/y/yozhao

  ## CLIC env ##
  
  ## setup latest GCC, CMake, ROOT, Geant4, Octave, Etc.
  
  source /cvmfs/clicbp.cern.ch/x86_64-centos7-gcc8-opt/setup.sh
  
  ## Variables for Geant4
  
  export G4_CMAKE_DIR=$G4INSTALL/lib64/Geant4-$(geant4-config --version)

## on lxclicbpk20 ##

elif [[ $(hostname) == "lxclicbpk20" ]]; then 

  #PATH=${PATH}:~/usr/bin

  source /home/yozhao/public/geant4.10.04.p02-install/bin/geant4.sh
  
  ## Variables for Geant4
  
  export G4_CMAKE_DIR=/home/yozhao/public/geant4.10.04.p02-install/lib64/Geant4-10.4.2

fi
