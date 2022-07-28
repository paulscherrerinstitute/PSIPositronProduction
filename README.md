![Logo of the project](Logo_Pcube.png)


# PSIPositronProduction

This repository collects the codes (Python, Matlab, ...) developed in relation to the PSI Positron Production (P-cubed) Project.
Most of the codes have been developed to

 * run simulations and
 * analyze and publish the related results.


## Getting started

### Prerequisites


#### For the analysis of simulation results

* Python3
* Jupyter lab (or Jupyter notebooks)

The latest versions can be easily installed with conda.
A conda environment file is available at `GIT_PSIPositronProduction/RepoSetup/Conda/JupyterNb.yml`:

```shell
$ cd GIT_PSIPositronProduction/RepoSetup/Conda
$ conda env create -f JupyterNb.yml
```

For more details, see:

https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file

#### To build Geant4 simulations

* Geant4, version 10.04.p03

See [How to build Geant4 on PSI RHEL7](Geant4OnPsiRhel7.md).


#### To edit codes with VSCode

The installation of VSCode on Linux or Windows is generally very simple.

This is not true if you are running a PSI RHEL7 Linux system. See [How to build VSCode from source on PSI RHEL7](VSCodeOnPsiRhel7.md) if you want (try) to keep going with VSCode on RHEL7...


### Clone the repository

**Warning**
In addition to a common installation of [Git](https://git-scm.com), you also need [Git Large File Storage](https://git-lfs.github.com) in order to download the large data files (usually compressed as `.tar.xz`) in the subfolder `GIT_PSIPositronProduction/Data`

Clone the repository with:

```shell
$ git clone git@github.com:paulscherrerinstitute/PSIPositronProduction.git GIT_PSIPositronProduction
```

or

```shell
$ git clone https://github.com/paulscherrerinstitute/PSIPositronProduction.git GIT_PSIPositronProduction
```

In this documentation, we indicate the root folder of the repository with `GIT_PSIPositronProduction`.


### Initial configuration

Some environment variables must be set. Follow the procedure for your operating system.


#### Linux

To temporarily set the environment variables:

  ```shell
  $ cd GIT_PSIPositronProduction/RepoSetup
  $ source Set_EnvVariables_Linux.sh
  ```

To permanently set the environment variables (recommended), add the following line to your `~/.bashrc` file:

  ```shell
  $ source PSIPositronProduction/RepoSetup/Set_EnvVariables_Linux.sh
  ```


#### Windows

Surf to `PSIPositronProduction/RepoSetup` and execute the batch script `Set_EnvVariables_Windows.bat` by double-clicking it.

**Warning**
This batch script does not check whether the desired paths are already present in the environment variables. Paths might be inserted multiple times in the same environment variable (should not be a problem). To avoid duplications, check the environment variables first:

  ```batch
  $ set
  ```

**Note**
It is a good practice to also do this after execution of the script. Important: open a NEW command window!



## Features


### Jupyter notebooks (.ipynb)

To run an analysis:

* Open the Jupyter notebook in JupyterLab (recommended) or JupyterNotebook.
* In the menu: Run > Run All Cells.
* Wait for the evaluation up to the last cell.


### FCCeeInjectorBeamApp

Use it directly through [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/paulscherrerinstitute/PSIPositronProduction/master?labpath=FCCeeInjectorBeamApp%2FFCCeeInjectorBeamApp.ipynb)


### Geant4 simulations

* [FCC-ee Target Tracking](FCCeeTargetTracking.md) (reference example)


### RF-Track simulations

(To be documented)


### Elegant simulations

(To be documented)


### ASTRA simulations

(to be documented)


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


## Deploying / Publishing

No active pipeline and no Git Pages at the moment.

