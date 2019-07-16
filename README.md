### micki

A modular, extensible, robust object-oriented microkinetic modeling package
written in Python.

### DEPENDENCIES:
 * lapack (MKL by default "-lmkl_rt"; can specify alternative LAPACK library via MICKI_LAPACK environmental variable)
 * ase
 * numpy
 * sympy
 * sundials (C library) - Tested with version 4.0. Compile Sundials with the following flags to cmake: -DFCMIX_ENABLE=ON -DCMAKE_C_FLAGS="-fPIC" -DLAPACK_ENABLE=ON -DSUNDIALS_INDEX_SIZE=32

