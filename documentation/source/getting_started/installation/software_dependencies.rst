.. ------------------------------------------------------------------------------
     (c) Crown copyright 2023 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------
.. _software dependencies:

Software dependencies of LFRic applications
===========================================

Compiler versions
-----------------

LFRic uses several features of modern Fortran including object
oriented features. Not all compilers have correctly implemented all
features. The following compilers are routinely tested at the Met
Office:

 * The Gnu compiler (versions 11.2 and 12.1)
 * The Intel ifort compiler (version 19.0)

Software Stack
--------------

To build and run typical LFRic applications, the following software
will be required. The numbers in parenthesis identify current versions
in use at the Met Office

Common software which may already be installed on HPC and research platforms

 * Python version 3 (3.7)
 * HDF5 (1.12.1)
 * NetCDF (4.8.1)
 * MPI (mpich 3.4.1)
 * blitz (1.0.2)

More specialist software essential for building and running LFRic applications.

 * PSyclone (2.5.0), a code generation library used by LFRic for
   generating portable performance code.
 * fparser (0.1.4), a Fortran parser used by PSyclone.
 * YAXT 0.10.0) which supports MPI data exchange in LFRic
   applications: https://swprojects.dkrz.de/redmine/projects/yaxt


Specialist software essential for fully-testing LFRic code developments.

 * PFUnit (3.3.3) A Fortran unit testing framework
 * stylist (0.4.1a): A code style-checker
