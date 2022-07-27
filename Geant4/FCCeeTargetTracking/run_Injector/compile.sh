#!/bin/bash

function do_compile(){

  module purge
  module load gcc/7.3.0
  module load geant4/10.5_multithreaded
  module load root/6.12.06
  source geant4.sh

  PD_RUN=$(pwd)
  
  cd ..
  
  PD=$(pwd)
  
  if [ $# -eq 0 ];then
    obj_name=Injector
  else
    obj_name=$1
  fi
  
  rm -rf ${obj_name}_build
  mkdir -pv ${obj_name}_build
  cd ${obj_name}_build
  
  NCore=$(nproc --all)
  
  cmake3 -DGeant4_DIR=${G4_CMAKE_DIR} ${PD}/${obj_name}
  
  make clean
  make -j ${NCore} ${obj_name}
  make install
  
  cd $PD_RUN

}

do_compile

