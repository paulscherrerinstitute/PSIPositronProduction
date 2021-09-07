# FCC-ee Target Tracking

This represents our reference example in Geant4.


## Building and executing the code

### On a personal Linux machine

Compile the code:

```shell
$ cd GIT_PSIPositronProduction/Geant4/FCCeeTargetTracking/run_Injector
$ ./compile.sh
```

Check the presence of the executable `../Injector_build/injector`.

Set the desired simulation parameters in the file `config.mac`.

Run the simulation:

```shell
$ ./run.sh
```

By default, the simulation output files are generated in the current directory (`FCCeeTargetTracking/run_Injector`).


### On Merlin6

For the general use of Merlin6, please refer to the official documentation:

https://lsm-hpce.gitpages.psi.ch/merlin6/introduction.html

First, login on a Merlin6 [login node](https://lsm-hpce.gitpages.psi.ch/merlin6/hardware-and-software.html#login-nodes), e.g. with NoMachine.

Load the necessary modules:

```shell
$ module load gcc/7.3.0
$ module load geant4/10.5_multithreaded
$ module load root/6.12.06
```

Compile the code like on a personal machine:

```shell
$ cd GIT_PSIPositronProduction/Geant4/FCCeeTargetTracking/run_Injector
$ ./compile.sh
```

To run the simulation, a job must be submitted to the queue through [Slurm](https://slurm.schedmd.com):

```shell
$ sbatch run_Merlin.sh
```

To check the progress of your job:

```shell
$ squeue
```


## Analyzing the results

The main output of the simulation is stored by default in the Root file `FCCeeTargetTracking/run_Injector/FCCeeTargetTracking.root`.

The results can be converted to the [Pcubed standard format](PcubedStandardFormat.md) e.g. with Python and the [PyROOT](https://root.cern/manual/python) package.

TODO: See `JupyterNb.ipynb`...
