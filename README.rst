*******
SkProgs
*******

Package containing a few programs that are useful in generating Slater-Koster
files for the DFTB-method.

**NOTE**: This packages comes with minimal documentation and with a currently
rather fragile user interface. It is considered to be neither stable nor
robust. Make sure, you check results as careful as possible. Use at your own
risk!


Installation
============

|build status|

Prerequisites
-------------

* Fortran 2003 compliant compiler

* CMake (>= 3.16)

* Python3 (>= 3.2)

* LAPACK/BLAS libraries (or compatible equivalents)

* libXC library with f03 interface (>=6.0.0)

* MpiFx (>=1.5, MPI-enabled build only)


Obtaining via Conda
-------------------

The preferred way of obtaining SkProgs is to install it via the conda package
management framework using `Miniconda
<https://docs.conda.io/en/latest/miniconda.html>`_ or `Anaconda
<https://www.anaconda.com/products/individual>`_. Make sure to add/enable the
``conda-forge`` channel in order to be able to access SkProgs::

  conda config --add channels conda-forge
  conda config --set channel_priority strict

We recommend to set up a dedicated conda environment and to use the
`mamba solver <https://mamba.readthedocs.io/>`_ ::

  conda create --name skprogs
  conda activate skprogs
  conda install conda-libmamba-solver
  conda config --set solver libmamba

to install the latest stable release of SkProgs (Fortran and Python
components)::

  mamba install skprogs skprogs-python


Building from source
--------------------

Follow the usual CMake build workflow:

* Configure the project, specify your compilers (e.g. ``gfortran``),
  the install location (i.e. path stored in ``YOUR_SKPROGS_INSTALL_FOLDER``,
  e.g. ``$HOME/opt/skprogs``) and the build directory (e.g. ``_build``)::

    FC=gfortran cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=YOUR_SKPROGS_INSTALL_FOLDER -B _build .

  An MPI enabled build is obtained by additionally setting ``-DWITH_MPI=1``
  (default: ``-DWITH_MPI=0``). At the moment only the two-center integration
  code ``sktwocnt`` is MPI parallelized and benefits from multiple processors.

  If libXC is installed in a non-standard location, you may need to specify
  either the ``CMAKE_PREFIX_PATH`` environment variable (if libXC was built with
  CMake) or the ``PKG_CONFIG_PATH`` environment variable (if libXC was built
  with autotools) in order to guide the library search::

    CMAKE_PREFIX_PATH=YOUR_LIBXC_INSTALL_FOLDER FC=gfortan cmake [...]

    PKG_CONFIG_PATH=FOLDER_WITH_LIBXC_PC_FILES FC=gfortran cmake [...]

* If the configuration was successful, build the code ::

    cmake --build _build -- -j

* After successful build, you should test the code by running ::

    pushd _build
    ctest -j
    popd

* If you want to test the MPI enabled binary with more than one MPI-process, you
  should set the ``TEST_MPI_PROCS`` variable in ``config.cmake`` accordingly,
  e.g.::

    set(TEST_MPI_PROCS "2" CACHE STRING "Nr. of processes used for testing")

  The ``TEST_MPI_PROCS`` cache variable can be updated or changed also after the
  compilation by invoking CMake with the appropriate ``-D`` option, e.g.::

    cmake -B _build -DTEST_MPI_PROCS=2 .
    pushd _build; ctest; popd

* If the tests were successful, install the package via ::

    cmake --install _build


Building libXC from source
--------------------------

Follow the usual CMake build workflow:

* Clone the official libXC repository and checkout the latest release tag, e.g.
  ``6.2.2``::

    git clone https://gitlab.com/libxc/libxc.git libxc
    cd libxc/
    git checkout 6.2.2

* Configure the project, specify your compilers (e.g. ``gfortran`` and ``gcc``),
  the install location (i.e. path stored in ``YOUR_LIBXC_INSTALL_FOLDER``, e.g.
  ``$HOME/opt/libxc``) and the build directory (e.g. ``_build``)::

      FC=gfortran CC=gcc cmake -DENABLE_FORTRAN=True -DCMAKE_INSTALL_PREFIX=YOUR_LIBXC_INSTALL_FOLDER -B _build .

* If the configuration was successful, build the code ::

    cmake --build _build -- -j

* After successful build, you should test the code by running ::

    pushd _build
    ctest -j
    popd

* If the tests were successful, install the package via ::

    cmake --install _build


Advanced build configuration
============================

Controlling the toolchain file selection
----------------------------------------

You can override the toolchain file, and select a different provided case,
passing the ``-DTOOLCHAIN`` option with the relevant name, e.g.::

  -DTOOLCHAIN=gnu

or ::

  -DTOOLCHAIN=intel

or by setting the toolchain name in the ``SKPROGS_TOOLCHAIN`` environment
variable. If you want to load an external toolchain file instead of one from the
source tree, you can specify the file path with the ``-DTOOLCHAIN_FILE`` option
::

  -DTOOLCHAIN_FILE=/path/to/myintel.cmake

or with the ``SKPROGS_TOOLCHAIN_FILE`` environment variable.

Similarly, you can also use an alternative build config file instead of
`config.cmake` in the source tree by specifying it with the
``-DBUILD_CONFIG_FILE`` option or by defining the ``SKPROGS_BUILD_CONFIG_FILE``
environment variable.


Generating SK-files
===================

The basic steps of generating the electronic part of the SK-tables are as
follows:

* If you have build SkProgs from source, initialize the necessary environment
  variables by sourceing the ``skprogs-activate.sh`` script (provided you have
  BASH or a compatible shell, otherwise inspect the script and set up the
  environment variables manually)::

    source <SKPROGS_INSTALL_FOLDER>/bin/skprogs-activate.sh

* Then create a file ``skdef.hsd`` containing the definitions for the elements
  and element pairs you wish to create. See the ``examples/`` folder for some
  examples.

* Run the ``skgen`` script to create the SK-tables. For example, in order to
  generate the electronic part of the SK-tables for C, H and O with dummy (zero)
  repulsives added, issue ::

    skgen -o slateratom -t sktwocnt sktable -d C,H,O C,H,O

  For an MPI enabled binary, make sure to prepend any required information to
  the two-center binary, e.g.::

    skgen -o slateratom -t "mpirun -np 2 sktwocnt" sktable -d C C |& tee output

  The SK-files will be created in the current folder. See the help (e.g. ``skgen
  -h``) for additional options.

Further documentation will be presented in a separate document later.


License
=======

SkProgs is released under the GNU Lesser General Public License.

You can redistribute it and/or modify it under the terms of the GNU Lesser
General Public License as published by the Free Software Foundation, either
version 3 of the License, or (at your option) any later version. See the files
`COPYING <COPYING>`_ and `COPYING.LESSER <COPYING.LESSER>`_ for the detailed
licensing conditions.


.. |build status| image:: https://img.shields.io/github/actions/workflow/status/dftbplus/skprogs/build.yml
    :target: https://github.com/dftbplus/skprogs/actions/
