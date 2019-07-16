### micki

A modular, extensible, robust object-oriented microkinetic modeling package
written in Python.

### DEPENDENCIES:
 * lapack (MKL by default "-lmkl_rt"; can specify alternative LAPACK library via MICKI_LAPACK environmental variable)
 * ase
 * numpy
 * sympy
 * sundials (C library) - Please use version 2.7.0 for best performance. Compatability with 4.X requires modifications to fortran.py (see comments). Compile with the following flags: -DFCMIX_ENABLE=ON -DCMAKE_C_FLAGS="-fPIC" -DLAPACK_ENABLE=ON
 

### model.py
Currently, the user must modify the compiler flags passed when compiling with f2py on lines 883-892 in model.py in order to use Micki. Current flags are what work for our system and are not guaranteed to work for abritrary systems.
