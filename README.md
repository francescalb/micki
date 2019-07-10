### micki

A modular, extensible, robust object-oriented microkinetic modeling package
written in Python. This software is not yet stable (but it's close!).

### DEPENDENCIES:
 * ase
 * numpy
 * sympy
 * sundials (C library) - has been tested successfully with versions 2.7.0 and 2.0.0. 

### model.py
Currently, the user must modify the compiler flags passed when compiling with f2py on lines 883-892 in model.py in order to use Micki. Current flags are what work for our system and are not guaranteed to work for abritrary systems. We will be working on eliminating this user modification to make it easier to use for future users.
