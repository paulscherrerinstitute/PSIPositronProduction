#!/bin/bash
#SBATCH --cluster=merlin6                 # Cluster name
#SBATCH --partition=hourly,daily,general  # Specify one or multiple partitions
#SBATCH --time=00:15:00                   # Strongly recommended
#SBATCH --hint=nomultithread              # Mandatory for multithreaded jobs
##SBATCH --exclusive                      # Uncomment if you need exclusive node usage
##SBATCH --ntasks-per-core=1              # Only mandatory for multithreaded single

## Advanced options example
##SBATCH --nodes=1                        # Uncomment and specify #nodes to use
##SBATCH --ntasks=44                      # Uncomment and specify #nodes to use
##SBATCH --ntasks-per-node=44             # Uncomment and specify #tasks per node
##SBATCH --cpus-per-task=44               # Uncomment and specify the number of cores per task
##SBATCH --output=<output_file>           # Generate custom output file
##SBATCH --error=<error_file>             # Generate custom error  file

module purge
module load gcc/7.3.0
module load geant4/10.5_multithreaded
module load root/6.12.06

source geant4.sh

config_file=config.mac
output_filename=FCCeeTargetTracking

## TODO: Understand how to use multiple cores on Merlin6
## ncore=$(nproc --all)
ncore=1
  
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
