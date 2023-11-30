TreElM
======

The actual sources for the treelm library are found in the tem-source repository.
It is included here in the `tem` subdirectory and amended with the environment to
build it along with aotus.

Use `git clone --recurse-submodules` when cloning this repository to fetch the
gathered subdirectories from the various repositories.

Prerequisite for building the library is an installed Python, Fortran compiler
and MPI library. For compilation you need to point `FC` to the appropiate MPI
compiler wrapper. (Usually `export FC=mpif90`).

The library can then be built with

```
bin/waf configure build
```

To install it, run:

```
bin/waf install
```

Run `bin/waf --help` to see all options.
