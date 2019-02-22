# LearnHatreeFock.jl
<img src="JuliaPsi.png" alt="JuliaPsi Logo" width="150"/>

## Learning Hatree-Fock By Doing

This package is focused on learning how to code a package for solving simple Hatree-Fock electronic structure calculations. 

The aim is to assist users of mainstream quantum chemistry and electronic structure software (e.g. NWChem, CP2K) with a more fundamental understanding of the steps associated with with all-electron wavefunction based approaches. This code will only address the Hatree-Fock approach (i.e. only electron exchange) using simple Gaussian basis-set (e.g. 6-31G). 

The code here is directly based on the implementation of Hatree-Fock in Matlab by [Prof. James Johns at the Univ. of Minnesota](http://www1.chem.umn.edu/groups/johns/index.html). This version is written in [Julia](https://julialang.org/), which will have familiar syntax as Matlab, but should have promising performance advantages and uses a liberal MIT license.

## Folder Hierarchy
  - /src - The primary source code for running HatreeFock.jl
  - /tests - Test modules for comparing to original code by Prof. James John.
  - /docs - Documentation.
  - /examples  - Example calculations and basis-sets.


