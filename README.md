# grasp-atomic-suite

All retained CSF consumers now use native descriptor Parquet + matching header TOML.
See [native input and build requirements](docs/implemented/design-0001-native-csf-descriptors.md).
The reader source is bundled in this repository and builds into `build/lib/libcsf_descriptor.a`;
no external rCSFs checkout or reader library is needed. Python is only used for developer tests.

[中文](README_ZH.md) | English

`grasp-atomic-suite` is a Fortran research and development repository for relativistic atomic-structure and atomic-property calculations. It combines original RMCDHF programs with `rhfs`, `ris4`, `rdensity`, `gj`, and their MPI support for hyperfine-structure, isotope-shift, and Landé `g_J` calculations.

This repository is not a minimal upstream GRASP mirror. It keeps the traditional GRASP numerical libraries and application layout, while adding CMake builds, MPI targets, an extracted `gj` program, regression data, and development notes for `g_J` and MPI optimization work.

Ordinary application directories now omit the old `90` suffix; the Landé-factor commands are `gj` and `gj_mpi`. See the [name mapping and rebuild instructions](docs/implemented/design-0005-program-naming.md).

## Source Layout

```text
.
├── configure.sh              # out-of-source CMake configuration helper
├── CMakeLists.txt            # top-level CMake build
├── src/
│   ├── apps/
│   │   ├── gj/               # standalone Landé g_J calculator
│   │   ├── rdensity/         # serial/MPI radial density and natural orbitals
│   │   ├── rhfs/             # relativistic hyperfine-structure program
│   │   ├── rmcdhf{,_mem,_mpi,_mem_mpi}/ # fresh GRASP copies with optional NP
│   │   └── ris4/             # relativistic isotope-shift program
│   └── lib/
│       ├── grasp_shared/     # shared state, capacities, and interfaces
│       ├── core/             # Shared CSF input, radial numerics, and utilities (libcore.a)
│       ├── davidson_eigensolver/ # Davidson eigensolver (libdavidson_eigensolver.a)
│       ├── breit_angular_support/ # Breit angular coefficients and external entry points
│       ├── spin_angular/        # angular algebra, Racah/tensor matrix elements
│       └── mpi/              # MPI file, path, and parallel helper routines
├── docs/                     # status index and classified documentation
│   ├── implemented/          # current implementation guides and code traces
│   ├── plans/                # experiments and validation still to run
│   ├── reviews/              # completed reviews and evidence
│   ├── reference/            # theory background drafts
│   └── archive/              # superseded plans and historical reports
├── data/                     # sample inputs and reference outputs for gj/RHFS checks
├── test/                     # orbital optimization and shared-library tests
├── scripts/                  # grid patcher and build cleanup helpers
├── bin/                      # installed executables
└── lib/                      # installed static libraries and Fortran modules
```

## Programs

| Program | Source directory | CMake targets | Purpose |
| --- | --- | --- | --- |
| Ordinary GRASP calculation chain | `src/apps/rangular*`, `rwfnestimate`, `rci*`, `rbiotransform*`, `rtransition*` | `rangular{,_mpi}`, `rwfnestimate`, `rci{,_mpi}`, `rbiotransform{,_mpi}`, `rtransition{,_phase,_mpi}` | Ten targets reuse the suite libraries and central capacities. See the [import record](docs/archive/2026-11-08-plan-0002-standard-grasp-integration.md). |
| `jj2lsj_omp` | `src/apps/jj2lsj` | `jj2lsj_omp` | Ordinary-CSF jj–LS conversion with OpenMP LS rows and ordered truncation; [build, runtime controls and validation](docs/implemented/design-0028-jj2lsj-openmp.md). CSFG support has been removed. |
| `rlevels` | `src/apps/rlevels` | Energy tables from RMCDHF `.mix` and RCI `.cimix` files with optional `jj2lsj` labels. [Usage and output format](docs/implemented/design-0025-rlevels.md). |
| `gj` | `src/apps/gj` | `gj`, `gj_mpi` | Standalone Landé `g_J` calculator extracted from the RHFS path. It reads `isodata`, the state's CSF descriptor pair, `name.mix/name.cimix`, and `name.wfn`, then writes `name.gj/name.cgj`. |
| `rhfs` | `src/apps/rhfs` | `rhfs`, `rhfs_mpi` | Relativistic hyperfine-structure program. It computes hyperfine constants and related matrix elements, writing `name.h/name.ch` and `name.hoffd/name.choffd`. |
| `rdensity` | `src/apps/rdensity` | `rdensity`, `rdensity_mpi` | Radial electron densities and statistically averaged natural orbitals; see [implementation and validation](docs/implemented/design-0026-rdensity-mpi.md). |
| `ris4` | `src/apps/ris4` | `ris4`, `ris4_mpi` | Relativistic isotope-shift program. It computes normal mass shift, specific mass shift, and field-shift electronic factors, writing `name.i/name.ci` and angular intermediate data. |
| `fical` | `src/apps/ris4` | `fical` | Standalone [line isotope-shift postprocessor](src/apps/ris4/README.md#fical) for pairs of RIS4 reports. |
| GRASP RMCDHF | `src/apps/rmcdhf`, `rmcdhf_mem`, `rmcdhf_mpi`, `rmcdhf_mem_mpi` | `rmcdhf`, `rmcdhf_mem`, `rmcdhf_mpi`, `rmcdhf_mem_mpi` | Copies from the workspace GRASP source, retaining its solver flow. Ordinary `rmcdhf_mpi` has no NP initialization; the other three optionally read `npol.inp`. See [NP interface and provenance](docs/implemented/design-0022-rmcdhf-nuclear-polarization.md). |

The rejected strict/diagnostic RMCDHF variants have been [removed](docs/implemented/adr-0002-retire-strict-and-diagnostic-rmcdhf.md), together with their dedicated tests and runtime entry points.

`./configure.sh` enables all retained non-CSFG applications by default: serial,
MPI and jj2lsj OpenMP, including `rwfnestimate`. MPI requires a compatible
Fortran MPI environment. Direct CMake configuration still requires
`-DGRASP_BUILD_SERIAL_APPS=ON` to include ordinary serial applications.
`OMP_NUM_THREADS` selects the OpenMP thread count; transformation-file modes in
`jj2lsj_omp` retain the serial backend.

## Build

Suite can build directly with GNU Fortran 16. The servers provide `gfortran`
16.2.1 (Red Hat 16.2.1-3). Load compatible MPI and BLAS/LAPACK first.
Run these commands from the Suite root and check that both the Fortran compiler
and Open MPI's `mpifort` report major version 16.

```fish
set -gx FC (command -v gfortran)
set -gx OMPI_FC $FC
$FC -dumpfullversion
mpifort -dumpfullversion
./configure.sh --build-dir build-gf16
cmake --build build-gf16 --parallel 4
```

Use `gfortran-16` if the server installs a versioned command. Choose a fresh
build directory when switching compilers. In the GraspKit workspace, put build
directories under `../temp/tests/`. For direct CMake configuration, pass
`-DCMAKE_Fortran_COMPILER=$FC`.

CMake enables only Fortran and uses it to find BLAS and LAPACK.
Cargo builds the bundled reader with Rust 1.99.0 or newer and uses a C
compiler for its Zstd dependency.
Cargo uses `--locked` and downloads missing dependencies by default.
Set `-DGRASP_CSF_READER_OFFLINE=ON` to build from the local cache.
`RCSFS_TEST_CLI` is only an offline test-fixture producer.
For CTest cases that generate descriptor fixtures, put a V3-capable `rcsfs`
on `PATH` or pass `-DRCSFS_TEST_CLI=/path/to/rcsfs` during configuration.

The default install prefix is `<build>/install`; use `--prefix PATH` or
`-DCMAKE_INSTALL_PREFIX=PATH` for another destination. A normal build does not
install or overwrite source `bin/`. Use `f2018-check` or `f2018-audit` to run
[strict acceptance on GNU Fortran 16](test/f2018/README.md).

The resulting executables are in `<build>/bin`, including `rwfnestimate`, the four
ordinary RMCDHF programs, MPI property programs, `rangular`, `rangular_mpi` and `jj2lsj_omp`.
CSFG stays off by default. Explicit `--no-serial`, `--no-mpi` or `--no-openmp`
options select a smaller build; `--serial` remains accepted for compatibility.

To build the four original RMCDHF programs and all property programs together:

```fish
cmake -S . -B build-all -DGRASP_ENABLE_MPI=ON -DGRASP_BUILD_SERIAL_APPS=ON \
  -DCMAKE_Fortran_COMPILER=$FC -DRCSFS_TEST_CLI=/path/to/rcsfs
cmake --build build-all --parallel 4
ctest --test-dir build-all --output-on-failure
cmake --install build-all
```

Debug build:

```fish
./configure.sh --debug
cmake --build build-debug --target install -j4
```

To disable MPI and build the serial applications, configure CMake directly:

```fish
mkdir build-serial
cmake -S . -B build-serial -DGRASP_ENABLE_MPI=OFF -DGRASP_BUILD_SERIAL_APPS=ON \
  -DCMAKE_Fortran_COMPILER=$FC -DRCSFS_TEST_CLI=/path/to/rcsfs
cmake --build build-serial --target install -j4
```

Rerunning `./configure.sh` updates an existing build from the same source tree
and compiler, including old caches where serial applications were disabled.
Defaults are applied again, so repeat opt-out options for a reduced build.
Use a new `--build-dir` when changing compilers or source trees. To remove build
artifacts explicitly, run:

```fish
./scripts/clean-build-artifacts.sh
```

## Dependencies And Configuration

The build requires a Fortran compiler, CMake 3.31 or newer, BLAS, and LAPACK. MPI targets additionally require an MPI Fortran environment such as `mpifort`.

CMake finds BLAS/LAPACK and requires GNU Fortran 15 or newer. The default strict
acceptance budget covers GNU Fortran 16.2 releases, including 16.2.1. The default
`f2018-release` profile uses O2; `f2018-check` and `f2018-audit` provide
standard, interface and runtime or optimized diagnostics. Legacy profiles,
implicit-SAVE flags and argument-mismatch allowances have been removed.
See the [build policy](docs/implemented/design-0002-fortran-build-profiles.md) and
[implemented compatibility and interface contracts](docs/implemented/design-0004-f2018-compatibility.md).
Current native input requirements and validation limits are documented in
[the descriptor guide](docs/implemented/design-0001-native-csf-descriptors.md).
Historical test counts do not certify the current source or physical accuracy.

## Inputs And Outputs

All interactive applications use the shared `Name of state:` prompt. Programs
that need two states label them as initial/final or, for FICAL, upper/lower.
MPI reads the answer on rank 0. Supported command-line forms continue to
accept a state directly. [Shared state input](src/lib/core/calculation_state_files.f90)

Run `rwfnestimate`, `rangular` and `rmcdhf` in the calculation directory and
enter the same name at each `Name of state:` prompt. For example, the state
`Zn_odd_AS1` uses `Zn_odd_AS1_desc.parquet` and `Zn_odd_AS1_header.toml`.

| Step | State files read | State files written |
| --- | --- | --- |
| `rwfnestimate` | `name_header.toml` | `name.initial.wfn` (initial radial orbitals) |
| `rangular` | `name_desc.parquet` and its header | `mcp.NN`, `name.rangular.log` |
| `rmcdhf` and `rmcdhf_mem` | Same descriptor pair, `name.initial.wfn`, `mcp.NN` | `name.wfn` (optimized orbitals), `name.mix` (OL/EOL mixing), `name.rmcdhf.summary`, `name.rmcdhf.log` |
| `rci` | Same descriptor pair and `name.wfn` | `name.cimix` (CI mixing), `name.csum`, `name.clog` |

MPI uses the same state name. Each `rangular_mpi` rank writes `mcp.mcpidx` in
its scratch directory. Run `rmcdhf_mpi` or `rmcdhf_mem_mpi` with the same basis
and rank layout. Use the current `rangular_mpi` to regenerate older MCP output.
The descriptor and radial files live in the permanent I/O directory, and the
angular and SCF logs live in the launch directory. `isodata` and applicable
`npol.inp` are shared nuclear inputs. See [MPI MCP formats and validation](docs/implemented/design-0020-rmcdhf-mcp-io.md).
RCI keeps its settings-first prompt order. Property programs such as `gj` use
up to 23 characters for the state name.

RMCDHF and RCI, including their MPI variants, require at least one valid ASF
root for every `J/parity` block in the CSF basis. Enter one line of root indices
per block. A blank line, a zero, or another invalid selection stops the run and
identifies the affected block. RMCDHF includes each selected root in the shared
orbital optimization and reads the EOL weight option when multiple roots are
selected. For legacy diagonal-only (E)AL MCP data, regenerate the angular
coefficients with the corresponding Suite `rangular` program.

MCP data belongs to the working directory. Finish the SCF run for a basis
before generating another basis there, or use separate work and MPI scratch
directories. States named `A` and `A.initial` also share the path `A.initial.wfn`
and need separate directories. State readers use `.wfn` for orbitals, `.mix` for
RMCDHF coefficients, and `.cimix` for RCI coefficients.

An empty state selects `rcsf_desc.parquet`, `rcsf_header.toml`, `rwfn.inp`,
`rwfn.out` and `rmix.out`. The legacy answer format starts with the y/n settings
response. For a single-letter state `Y` or `N` (either case), the
matching header in the input directory takes precedence. Fixed-name runs can
still be archived with [save_rmcdhf](scripts/README.md#保存-rmcdhf-计算结果).

Common input files:

- `isodata`: nuclear and isotope data.
- `name_desc.parquet` and `name_header.toml`: CSF descriptors and matching header.
- `name.wfn`: radial orbital wavefunctions in the [radial wavefunction stream](#radial-wavefunction-stream).
- `name.mix`: MCDHF/non-CI mixing-coefficient file.
- `name.cimix`: CI mixing-coefficient file.

Common output files:

- `gj`: `name.gj` or `name.cgj`.
- `rhfs`: `name.h/name.ch` and `name.hoffd/name.choffd`.
- `rdensity`: `name.nw` and `name.d/name.cd`.
- `ris4`: `name.i/name.ci`, and possibly angular intermediate files such as `name.IOB` and `name.ITB`.

`gj_mpi` accepts a state name with `--nonci` or `--ci`. Run it from the
calculation directory and set `suite_bin` to the absolute path of the build's
`bin` directory. CI mode reads `name.cimix`.

```fish
mpiexec -n 2 "$suite_bin/gj_mpi" name --nonci
mpiexec -n 2 "$suite_bin/gj_mpi" name --ci
```

Without arguments, it falls back to the interactive input flow.

### Radial wavefunction stream

The `rwfnestimate` and RMCDHF wavefunction outputs use the `RAD--WFN` raw binary stream.
Integers are little-endian signed 32-bit values and floating-point values are little-endian IEEE binary64.
Offsets below are zero-based bytes. Suite readers also accept earlier `GRWFRAD1` streams and sequential Fortran `G92RWF` files.

| Byte offset | Data |
| --- | --- |
| 0 | Eight ASCII bytes `RAD--WFN` |
| 8 | Integers `K` (orbital count) and `L` (largest orbital point count) |
| 16 | `K` consecutive 32-byte orbital entries, each containing `n`, `kappa`, `point_count`, a reserved zero integer, `energy`, and `origin_parameter` in that order |
| `16 + 32*K` | `float64` matrix with C/NumPy shape `(2*K+1, L)` and Fortran shape `(L, 2*K+1)` |

The first matrix row contains the shared radial coordinates in bohr. The following rows contain `P` then `Q`
for each orbital in table order. Only the first `point_count` samples of each orbital are active and the rest are zero.
Orbital energies use the Suite binding-energy convention in hartree.

```python
import struct

import numpy as np

filename = "state.wfn"
with open(filename, "rb") as stream:
    magic, count, longest = struct.unpack("<8sii", stream.read(16))
    if magic not in (b"RAD--WFN", b"GRWFRAD1"):
        raise ValueError("Unknown wavefunction format")

orbitals = np.memmap(
    filename, mode="r", offset=16, shape=(count,),
    dtype=[
        ("n", "<i4"), ("kappa", "<i4"), ("point_count", "<i4"),
        ("reserved", "<i4"), ("energy", "<f8"), ("origin_parameter", "<f8"),
    ],
)
radial = np.memmap(
    filename, mode="r", dtype="<f8", offset=16 + 32 * count,
    shape=(2 * count + 1, longest),
)
index = 0
points = int(orbitals[index]["point_count"])
r = radial[0, :points]
p = radial[2 * index + 1, :points]
q = radial[2 * index + 2, :points]
```

### Indexed ASF mixing stream

RMCDHF `.mix` and RCI `.cimix` use the same `CSF--MIX` format. The stream uses little-endian signed
64-bit integers and IEEE binary64 floating-point values. Offsets and table row numbers are zero-based.
The descriptor SHA-256 is 64 lowercase ASCII hex digits.

| Region | Encoding |
| --- | --- |
| File header | Eight ASCII bytes `CSF--MIX`, five integers (electrons, orbitals, blocks, total CSFs, total ASF roots), then the descriptor SHA-256 |
| Block directory | One 48-byte entry per descriptor block. The six integers are physical `2J`, parity `+1/-1`, CSF count, root count, first row in the level table, and byte offset of this block's coefficients |
| Level table | One 16-byte entry per root, with its original one-based root number inside the block and absolute energy in hartree |
| Coefficients | Signed binary64 values in descriptor order. All CSFs for one root are contiguous. A block has NumPy/C shape `(n_roots, n_csfs)` and Fortran shape `(n_csfs, n_roots)` |

The block directory starts at byte 112 and the level table starts at `112 + 48*n_blocks`. The position
of a root's first coefficient is the block offset plus `8*n_csfs*root_row`, where `root_row` is its
zero-based position in the block's level table. Readers can use the file tables without reading the
coefficient payload. `GRMIXTMP` marks an incomplete write. The format readers accept earlier
`GRMIX001` streams and sequential Fortran `G92MIX` files.

```python
from struct import unpack

import numpy as np

filename = "state.cimix"
block_index = 0
with open(filename, "rb") as stream:
    magic, electrons, orbitals, blocks, csfs, levels, basis_sha = unpack(
        "<8s5q64s", stream.read(112)
    )
    if magic not in (b"CSF--MIX", b"GRMIX001") or not 0 <= block_index < blocks:
        raise ValueError("Unknown mixing format or block index")
    stream.seek(112 + 48 * block_index)
    two_j, parity, n_csfs, n_roots, first_level, coefficients_offset = unpack(
        "<6q", stream.read(48)
    )

level_table = np.memmap(
    filename, mode="r", offset=112 + 48 * blocks + 16 * first_level,
    dtype=[("root", "<i8"), ("energy_hartree", "<f8")], shape=(n_roots,),
)
coefficients = np.memmap(
    filename, mode="r", offset=coefficients_offset,
    dtype="<f8", shape=(n_roots, n_csfs),
)
first_root = level_table[0]
first_root_coefficients = coefficients[0]
```

## Verification Data

The `data/` directory contains a `test` data set:

```text
data/isodata
data/test.c
data/test.m
data/test.w
data/test.gj
```

To reproduce `data/test.gj`, run the following from the Suite root with a
V3-capable `rcsfs` on `PATH`. Use your actual build directory in `suite_bin`.
`rcsfs convert` creates both the Parquet descriptor and its matching TOML
header. Copy the sample orbitals and mixing coefficients to `test.wfn` and `test.mix`.

```fish
set suite_root (pwd)
set suite_bin (realpath build-gf16/bin)
mkdir -p ../temp/tests
set gj_work (mktemp -d ../temp/tests/gj-example.XXXXXXXX)
cp data/isodata data/test.c "$gj_work/"
cp data/test.m "$gj_work/test.mix"
cp data/test.w "$gj_work/test.wfn"
cd "$gj_work"
rcsfs convert test.c test_desc.parquet
printf '%s\n' y test n | "$suite_bin/gj"
diff -u "$suite_root/data/test.gj" test.gj
cd "$suite_root"
```

CTest covers descriptor input, calculation workflows, shared-library routines,
sparse MPI buffers, and the grid patcher.

```fish
ctest --test-dir build-all --output-on-failure
```

Orbital regression runners live in [test/rmcdhf](test/rmcdhf/README.md). Large calculation inputs/results remain outside the repository. The retired orbopt runners and SLURM templates are no longer provided; see [retirement details](docs/implemented/adr-0001-retire-rmcdhf-orbopt.md).

## Radial Grid And Migration

Shared limits live in [`cal_constants.f90`](src/lib/grasp_shared/cal_constants.f90)
and are re-exported by [`grasp_parameters.f90`](src/lib/grasp_shared/grasp_parameters.f90).
See [shared library naming](docs/implemented/design-0006-shared-library-naming.md) for the library and module mapping.
The Suite defaults are `max_radial_points=590`, `radial_workspace_capacity=600`, `max_orbitals=127`,
`max_orbitals_minus_one=126`, `max_orbitals_minus_two=125`, and `orbital_index_radix=215`.
Finite nuclei use `radial_point_count=max_radial_points`. Programs with a point-nucleus branch start at
`min(220,max_radial_points)`.
Input orbitals are interpolated onto the current calculation grid. Upstream RMCDHF initializes `ACCY=H**6` before optional
interactive grid changes; explicitly supply ACCY if you change H.
Copy [scripts/config.example.toml](scripts/config.example.toml) to `scripts/config.toml`,
edit the local settings, and run [scripts/patch_grasp_grid.py](scripts/README.md)
without arguments. The patcher locates `src/` alongside its own `scripts/` directory.
Both CLI `--nnnp` and TOML `nnnp` set the
compile-time radial-grid capacity `max_radial_points`.
Local `config.toml` is ignored by Git; omitted grid keys preserve existing values.
Preview is the default. Keep a separate source copy for each parameter set and
build it manually using this repository's CMake workflow. Changes require a
full rebuild and numerical convergence checks.

See [shared parameter guide](docs/implemented/design-0003-shared-grid-parameters.md) for the parameter table,
initialization order, and validation. The suite layout edits the central capacity configuration
and the fifteen individual grid entry points as needed. Explicit changes to `working_tolerance`
also cover the serial and MPI RCI restart paths.
For candidate values, their rationale, test datasets, and convergence criteria,
see [GRASP grid parameter tuning](docs/plans/plan-0001-radial-grid-validation.md) (Chinese).

The [RMCDHF source snapshot](docs/reference/rmcdhf_upstream_snapshot.json) records the source provenance. The [RMCDHF flow guide](docs/implemented/design-0019-rmcdhf-flow-and-state.md) describes the retained solver.

## Development Notes

The [documentation index and status audit](docs/README.md) covers every document,
its implementation evidence, superseded claims, and remaining validation work.

- [Implemented guides](docs/implemented/design-0003-shared-grid-parameters.md): shared parameters and the [RHFS g_J code trace](docs/implemented/design-0027-rhfs-lande-g.md).
- [Experiment plan](docs/plans/plan-0001-radial-grid-validation.md): grid configurations and physical convergence checks still to run.
- [Completed grid review](docs/reviews/evidence-0001-radial-grid-methodology-review.md): source findings and diagnostic evidence.
- [Theory reference](docs/reference/design-0032-lande-g-theory.md): Landé-factor background draft; operator conventions need independent verification.
- [Historical archive](docs/README.md#历史背景): gj blueprint, orbital-optimization results and RHFS/RIS optimization reports.

## Contributor Notes

- Preserve the local Fortran style of the file being edited; avoid broad formatting-only rewrites.
- Name application source files and Fortran symbols in lowercase `snake_case` by responsibility. Use the implementation name for migrated files whose modules retain a `*_interface` suffix. In `src/lib/grasp_shared/`, use one module per matching file. [Repository guidelines](AGENTS.md#coding-style--naming-conventions)
- Numerical changes should include input data, output differences, or regression notes.
- MPI changes should be checked against the corresponding serial path for input, output, and numerical consistency.

## License

This repository is distributed under the [MIT License](LICENSE).

## GRASPG / CSFG extension

The CSFG extension source, targets and dedicated tests were removed on 2026-11-16.
`GRASP_BUILD_CSFG=ON` is rejected. `--no-csfg` and `GRASP_BUILD_CSFG=OFF` remain
accepted for compatibility with existing ordinary-suite build commands.
The ordinary shared libraries remain in use. See the
[GRASPG source manifest](docs/reference/graspg/source-manifest.json) for provenance.

### Attribution and citation

The CSFG implementation is derived from [compas/graspg](https://github.com/compas/graspg).
This repository previously integrated and adapted that upstream code; the GRASPG implementation
and CSFG methods are credited to their original authors. Research using these
programs should cite the GRASPG software paper, the CSFG method paper, and the
underlying GRASP2018 package:

1. Ran Si, Yanting Li, Kai Wang, Chongyang Chen, Gediminas Gaigalas, Michel Godefroid, and Per Jönsson. “Graspg – An extension to Grasp2018 based on configuration state function generators.” *Computer Physics Communications* **312**, 109604 (2025). [DOI: 10.1016/j.cpc.2025.109604](https://doi.org/10.1016/j.cpc.2025.109604).
2. Yanting Li et al. “Reducing the computational load – atomic multiconfiguration calculations based on configuration state function generators.” *Computer Physics Communications* **283**, 108562 (2023). [DOI: 10.1016/j.cpc.2022.108562](https://doi.org/10.1016/j.cpc.2022.108562).
3. C. Froese Fischer, G. Gaigalas, P. Jönsson, and J. Bieroń. “GRASP2018—A Fortran 95 version of the General Relativistic Atomic Structure Package.” *Computer Physics Communications* **237**, 184–187 (2019). [DOI: 10.1016/j.cpc.2018.10.032](https://doi.org/10.1016/j.cpc.2018.10.032).
