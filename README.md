![Logo of the project](Logo_Pcube.png)


# PSIPositronProduction

This repository collects the codes (Python, Matlab, ...) developed in relation to the PSI Positron Production (Pcube) Project.
An important part of the codes has been developed to run simulations and to analyze and publish the related results.


## Getting started

### Prerequisites

#### To build Geant4 simulations

* Geant4, version 10.04.p03

See [How to build Geant4 on PSI RHEL7](Geant4OnPsiRhel7.md).


#### For the analysis of Elegant simulations

* Python3
* Jupyter notebooks

The latest versions can be easily installed with conda.
A conda environment file is available at `GIT_PSIPositronProduction/RepoSetup/Conda/JupyterNb.yml`:

```shell
$ cd GIT_PSIPositronProduction/RepoSetup/Conda
$ conda env create -f JupyterNb.yml
```

For more details, see:

https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file


#### VSCode on PSI RHEL7

See [How to build VSCode from source on PSI RHEL7](VSCodeOnPsiRhel7.md) if you want to keep going with VSCode on RHEL7.


### Clone the repository

```shell
$ git clone git@gitlab.psi.ch:schaer_m/psipositronproduction.git GIT_PSIPositronProduction
```

We indicate the root folder of the repository with GIT_PSIPositronProduction.


### Initial configuration

Set the environment paths to declare the location of the different Python modules:

* Linux:

  ```shell
  $ cd GIT_PSIPositronProduction/RepoSetup
  $ source Set_Pythonpath.sh
  ```

* Windows: not implemented yet.
* Mac: not implemented yet.


## Features

### Geant4 simulations

* [FCC-ee Target Tracking](FCCeeTargetTracking.md) (reference example)


### Elegant simulations

(To be documented)


## Code documentation

(To-do)


## Unit testing

(To-do)


## Contributing

Create your own branch and start contributing:

```shell
$ git checkout -b MyBranch
```

Please contact mattia.schaer@psi.ch to discuss the details.


### Deploying / Publishing

No active pipeline and no Git Pages at the moment.

