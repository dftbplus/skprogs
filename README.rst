*******
SkProgs
*******

Package containing a few programs enabling to generate Slater-Koster files for
the DFTB-method.

**NOTE**: This packages comes with minimal documentation and with a currently
rather fragile user interface. It is considered to be neither stable nor
robust. Make sure, you check results as careful as possible. Use at your own
risk!


Installing
==========

Prerequisites
-------------

* Fortran 2003 compiler

* CMake (>= 3.16)

* Python3

* LibXC library with f90 interface (tested with version 4.3.4, version 5.x does
  not work due to inteface changes in LibXC)

  
Building the code
-----------------

Follow the usual CMake build workflow:

* Configure the project, specify your compiler (e.g. ``gfortran``), the install
  location (e.g. ``$HOME/opt/skprogs``) and the build directory
  (e.g. ``_build``)::

    FC=gfortran cmake -DCMAKE_INSTALL_PREFIX=$HOME/opt/skprogs -B _build .

  If libXC is installed in a non-standard location, you may need to specify
  either the ``CMAKE_PREFIX_PATH`` environment variable (if libXC was built with
  CMake) or the ``PKG_CONFIG_PATH`` environment variable (if libXC was built
  with autotools) in order to guide the library search::

    CMAKE_PREFIX_PATH=YOUR_LIBXC_INSTALL_FOLDER FC=gfortan cmake [...]
    
    PKG_CONFIG_PATH=FOLDER_WITH_LIBXC_PC_FILES FC=gfortran cmake [...]

* If the configuration was successful, build the code ::

    cmake --build _build -- -j

* If the build was successful, install the code ::

    cmake --install _build


Generating SK-files
===================

The basic steps of generating the electronic part of the SK-tables are as
follows:

* Initialize the necessary environment variables by sourceing the
  ``skprogs-activate.sh`` script (provided you have BASH or compatible shell,
  otherwise inspect the script and set up the environment variables manually)::

    source <SKPROGS_INSTALL_FOLDER>/bin/skprogs-activate.sh

* Then create a file ``skdef.hsd`` containing the definitions for the elements
  and element pairs you wish to create. See the ``examples/`` folder for some
  examples.

* Run the ``skgen`` script to create the SK-tables. For example, in order to
  generate the electronic part of the SK-tables for C, H and O with dummy (zero)
  repulsives added, issue ::

    skgen -o slateratom -t sktwocnt sktable -d C,H,O C,H,O

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
