## RIS4: A program for relativistic isotope shift calculations - official development version
Note: In modern Fortran and compatible with GRASP2018. First version provided by Jörgen Ekman (Malmö University) on June 7th, 2020.

### Abstract
Spectral lines from different isotopes display a small separation in energy, commonly referred to as the line isotope shift. The program ris 4 (Relativistic Isotope Shift) calculates normal and specific mass shift parameters as well as field shift electronic factors from relativistic multiconfiguration Dirac–Hartree–Fock wave functions. These quantities, together with available nuclear data, determine isotope-dependent energy shifts. Using a reformulation of the field shift, it is possible to study, in a model-independent way, the atomic energy shifts arising from changes in nuclear charge distributions, e.g. deformations.

Original publication: [Ekman et al. Comput. Phys. Comm. 235 (2019) 433](https://doi.org/10.1016/j.cpc.2018.08.017)

### FICAL

`fical` reads the two states' RIS4 `.i` or `.ci` reports and an `isodata` file in
the working directory. It prompts for the state names and writes
`<state1>.<state2>.fi` or `<state1>.<state2>.cfi`. It uses the unit on each
converted NMS and SMS row and reports differences in the requested GHz, MHz,
or meV units. A missing or unsupported mass-shift unit is an input error.
The linear-dependence discrepancy is reported as `NaN` when its normalization denominator is zero.

The suite builds and installs `fical` as a separate executable. You can
also compile it directly from the repository root with GNU Fortran 15.

```fish
/usr/bin/gfortran-15 -std=f2018 -pedantic-errors -Werror=implicit-interface -Werror=implicit-procedure src/apps/ris4/fical.f90 -o fical
```
