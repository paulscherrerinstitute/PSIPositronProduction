# How to build Geant4 on PSI RHEL7

Currently (11.08.2021), the official Linux system at PSI is Red Hat Enterprise Linux Server 7.9 (RHEL7).
In this page we provide specific, step-by-step instructions to install Geant4 on the mentioned system.

## On this page

[[_TOC_]]


## Reference page

The reference page for Geant4 is:

https://geant4.web.cern.ch

Instructions to build and run from source are at:

https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/installguide.html#buildandinstall

Note that Geant 4 can also be installed via several package managers, e.g. `conda`:

https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html


## Specific prerequisites for RHEL7

### Recent GCC toolchain

Install a recent GCC toolchain (e.g. `devtoolset-10`) and activate it:

```shell
$ sudo yum install devtoolset-10
$ scl enable devtoolset-10 bash
```

For more details, see:

https://access.redhat.com/documentation/en-us/red_hat_developer_toolset/10/html/user_guide

or

https://developers.redhat.com/products/developertoolset/hello-world


### Cmake3

Install `cmake3`:

```shell
$ sudo yum install cmake3
```

This package provides the `ccmake3` command (note the double "c") that we use to configure the build.


## Install Geant4

### Download version 10.4.3

We have chosen this version in order to be able to compile the example 'Geant4/FcceeTarget_StartingExample' as provided by CERN, without any modification.

Downloaded from
https://geant4.web.cern.ch/support/download_archive
(surf through the pages at the bottom to find the desired version).

Extract the archive, jump into the extracted folder and create the build directory:

```shell
$ tar -xf geant4.10.04.p03.tar.gz
$ cd geant4.10.04.p03
$ mkdir build
$ cd build
```


### Configure the build


First of all, do not forget to activate devtoolset!

From the directory `geant4.10.04.p03/build`:

```shell
$ ccmake3 ..
```

Select the desired build options and generate (key "g", available after having configured twice with key "c").
In the following you find the choice on my machine MPC1452:

```shell
CMAKE_BUILD_TYPE                 Release
CMAKE_INSTALL_PREFIX             /usr/local/geant4-10.4
GEANT4_BUILD_MULTITHREADED       ON
GEANT4_INSTALL_DATA              ON
GEANT4_INSTALL_DATADIR
GEANT4_USE_G3TOG4                ON
GEANT4_USE_GDML                  ON
GEANT4_USE_INVENTOR              OFF
GEANT4_USE_OPENGL_X11            ON
GEANT4_USE_QT                    ON
GEANT4_USE_RAYTRACER_X11         ON
GEANT4_USE_SYSTEM_CLHEP          OFF
GEANT4_USE_SYSTEM_EXPAT          ON
GEANT4_USE_SYSTEM_ZLIB           OFF
GEANT4_USE_XM                    ON
Qt5Core_DIR                      /usr/lib64/cmake/Qt5Core
Qt5Gui_DIR                       /usr/lib64/cmake/Qt5Gui
Qt5OpenGL_DIR                    /usr/lib64/cmake/Qt5OpenGL
Qt5PrintSupport_DIR              /usr/lib64/cmake/Qt5PrintSupport
Qt5Widgets_DIR                   /usr/lib64/cmake/Qt5Widgets
```

Depending on your choices, you will need to install additional packages with `yum` on your system.

If you intend to install several versoins in your system, it is recommended to change `CMAKE_INSTALL_PREFIX` from the default `/usr/local` to something more specific, e.g. `/usr/local/geant4-10.4` (like above).


### Compile and install

```shell
$ make -j 8
$ sudo make install
```

With the option `-j 8` you compile in parallel e.g. on 8 cores.


### Activate installation

```shell
$ source /usr/local/geant4-10.4/bin/geant4.sh
```


### Build an example

Build the Geant4 official example B1:

```shell
$ cd geant4.10.04.p03/examples/basic/B1
$ mkdir build
$ cd build
$ ccmake3 ..
(Configure and generate)
$ make
```

Execute it:

```shell
$ ./exampleB1
```

