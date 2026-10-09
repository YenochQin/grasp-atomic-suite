# RDENSITY serial and MPI migration

Status: implemented, 2026-10-09. Migration base: suite `47fd829`; upstream
`grasp/src/appl/rdensity` at `488e875d`. The migration described here builds on
those revisions. Local regression is separate from production model convergence
and large-CSF scaling validation.

## Build and run

The normal suite build includes `rdensity_mpi` when MPI is available. To also
build the serial reference and enable their comparison tests:

```sh
cmake -S . -B build-all -DGRASP_BUILD_SERIAL_APPS=ON
cmake --build build-all --parallel 1
ctest --test-dir build-all -R rdensity --output-on-failure
cmake --install build-all
```

The full local Debug build passed with one compilation job. Building all existing
serial/MPI applications together with four jobs exposed a pre-existing RHFS/g_J
module-file race (`gethfd_i.mod0`). The two new RDENSITY targets use separate
module directories and were built successfully with four compilation jobs.

Run in a directory containing `isodata`, `name.c`, `name.w`, and `name.m` or
`name.cm`. All ranks must see the same read-only input files.

```sh
mpiexec -n 4 /path/to/grasp-atomic-suite/bin/rdensity_mpi name --ci --order=2
/path/to/grasp-atomic-suite/bin/rdensity name --nonci --order=1
```

The default mixing suffix is `.m`; `--ci` selects `.cm`. The default ordering is
2. State names follow the legacy limit of 23 characters, with no spaces.
Without arguments, both targets accept the upstream default-input sequence:

```text
y
name
y
2
```

These answers select defaults, the state name, CI coefficients, and ordering.
Multi-rank runs support default settings, following the RHFS/RIS entry points.
Serial runs retain the upstream debug and custom radial-grid prompts. Both
targets use `radial_grid_defaults`: finite nuclei default to `N=NNNP=2990`,
not upstream RDENSITY's 590-point capacity. Use the
[common-parameter contract](common_parameters.md) to configure other grids.

Only rank zero writes `name.nw` and `name.d`/`name.cd`. Existing density files
cause an error before calculation; `.nw` is replaced, as in upstream. Occupations
and rotation coefficients are printed in output orbital order. The unformatted
`.nw` file keeps the `G92RWF` record structure and original orbital energy fields;
these fields are not natural-orbital energies. CI coefficients are not rotated
or written: perform a subsequent RCI calculation in the new orbital basis when
new mixing coefficients are needed.

## Parallel and numerical implementation

- [rdensity.f90](../../src/appl/rdensity/rdensity.f90) is shared by both targets.
  Only rank zero reads interactive options, then broadcasts them. Each rank
  independently loads the read-only CSFs, nuclear data, wavefunctions and mixing
  coefficients using the existing suite libraries. No library was copied.
- [natorbnew.f90](../../src/appl/rdensity/natorbnew.f90) assigns columns cyclically,
  `IC=rank+1,NCF,workers`, as in RHFS/RIS. Each process has its own mutable angular
  workspace. Exact J/parity and one-body occupation selection rules reject zero
  contributions before `ONESCALAR`; the coefficient cutoff remains `1D-10`.
- Each rank accumulates `rho(NW,NW,NVEC)`. `GDRSUMMPI_ROOT` reduces one state at a
  time. Idle ranks participate with zero matrices. Root computes radial densities
  in the original orbital basis, forms the `(2J+1)` statistical average over all
  input levels, and diagonalizes each actual kappa block with LAPACK `DSYEV`.
- The contraction uses diagonal terms plus twice each off-diagonal term. This is
  algebraically equivalent to upstream accumulation over CSF pairs; its floating
  point summation order differs. The output retains upstream's convention
  `D(r)=r^2*rho(r)` (its printed `rho` has no extra `1/(4*pi)` factor).
- The radial contraction happens once after reduction. The original
  `DINT1VEC(NNNW,NNNW,NNNP)` allocation and per-CSF-pair grid loop are removed.
  At the suite capacities that cache alone would occupy about 368 MiB **per
  rank**. Density-matrix storage now scales as `8*NW*NW*NVEC` bytes per rank;
  CSF and CI arrays are still replicated. The CSF pair traversal remains
  quadratic, and I/O, eigensystems and output remain serial on root.
- [natural_orbitals.f90](../../src/appl/rdensity/natural_orbitals.f90) shares the
  eigensolver and radial rotation between serial/MPI. Arrays use actual block
  sizes; nonconsecutive principal quantum numbers and more than 20 orbitals in
  a kappa block are supported. LAPACK status is checked. Ordering 1 requires
  unique dominant components. Ordering 2 reserves pure orbitals, then assigns
  remaining eigenvectors in decreasing occupation order without overwriting
  slots. P, Q and PZ phases are changed together after all PZ values are rotated.
  These correct the upstream fixed-size, ordering and PZ-index limitations.
- Local allocatable workspaces use recursive procedures so they have automatic
  lifetime despite the suite's `-fno-automatic` compiler flag. No OpenMP regions
  surround the stateful angular library.

Input dimension, truncated mixing-record, finite-coefficient and vector-norm
checks fail through a serial error or `MPI_Abort`. The older shared CSF/radial
readers retain their existing error handling. The migration does not establish
bitwise MPI reproducibility; sums may differ in their final rounding bits.
Degenerate natural-orbital subspaces are not individually unique.

## Validation

Local results: all 16 Debug CTest cases passed with array bounds checks enabled;
all five RDENSITY tests passed in Release, including 1/2/4-rank MPI; the separate
MPI-disabled Release build passed both RDENSITY tests. These are local regression
results, not a cluster speedup measurement.

[CTest registration](../../test/CMakeLists.txt) and
[regression driver](../../test/rdensity/regression.py) cover:

- The production eigensolver with a 23-orbital block, mixed/pure ordering,
  orthogonality, occupation eigenvalues, repeated calls, dominance ambiguity,
  and consistent P/Q/PZ phase changes.
- Generated one-electron GRASP files with analytically known diagonal and
  off-diagonal densities; two different J/parity blocks verify the 2:4 weights.
  Generated inputs are synthetic test fixtures, not physical reference models.
- Reconstruction of the weighted radial density from the written natural
  orbitals and their occupations, including a gap in principal quantum numbers.
- Serial versus 1-, 2- and 4-rank output comparisons for the synthetic cases
  and the existing Ni fixture (1042 CSFs, five J blocks, nine levels).
  Four ranks exercise ranks with no assigned columns in the small fixtures.
- Interactive ordering 1, CI/non-CI inputs, help, invalid arguments, missing
  files and truncated mixing data. All generated files use temporary directories.

These local tests check implementation consistency and mathematical invariants.
Large active-space performance, cluster scaling and scientific convergence of
production systems still require representative calculations.

An additional local comparison used the original upstream executable and the
synthetic two-J fixture on the same explicitly selected 590-point grid. P/Q
values agreed up to orbital phase within `8.5e-15`; printed densities agreed
except for subnormal rounding (maximum common-row absolute difference below
`2.1e-318`). The reordered contraction emitted one extra positive tail row of
`2.207e-320`. Exact text identity of positive-only density files is therefore
not a suitable regression criterion.
