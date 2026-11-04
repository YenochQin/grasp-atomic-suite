# jj2lsj_omp conversion regression

The fixtures are generated, small CSF bases, not atomic calculation results.
`check.py` creates native four-byte-record-marker GRASP mixing files with
four-byte integers and eight-byte reals using the Python standard library.
No solver, generator binary or external Python package is required at test time.

The `.c` fixtures were generated and expanded with the suite CSFG tools at base
revision `a4661441492e38beb99293182cb9257326d2dcef`:

| Fixture | Reference configurations | Maximum 2J | CSFs |
| --- | --- | --- | --- |
| even | 2p², 3p², 3d², 4d², 4f², 5f², 6h², 7h² | 12 | 82 |
| odd | 2p 3d, 4f 5g, 6h 7i | 12 | 60 |
| f3 | 4f³ | 9 | 29 |

Generator inputs use core option `0`, each listed configuration with `(N,i)`,
the listed orbital set, minimum 2J=0 (even/odd) or 1 (f3), zero excitations,
and `n` for another generation. The expanded ordinary CSFs are checked in.

Mixing vectors are either identity columns (complete transformation) or
normalized vectors concentrated in the first two jj-CSFs with small signed
components elsewhere. The latter ensure the other configurations survive
screening while the LS cumulative norm reaches its stopping threshold early.
They need not be physical eigenstates or mutually orthogonal in truncated cases.

`coefficients.f90` uses the production input reader, LS generator and conversion
dispatcher and dumps every signed double, including the truncated zero tail.
The script compares scalar results with 1/2/4/4-thread launches, using separate
processes to exercise initially cold angular caches. Full identity columns also
check transformation orthonormality to `2e-12`; scalar/parallel comparisons have
no tolerance and must be bitwise identical. The even fixture has blocks larger
than both two- and four-thread batches and exercises partially filled last batches.

The tests cover default settings, full transformation, 25% maximum omitted norm,
selected levels in reverse order with skipped blocks, sequential `.lsj.T`
generation/replay, and unique labels. Default `.expected.lbl` files were generated
with the original workspace GRASP executable, whose source snapshot is identified
in [the provenance manifest](../../docs/reference/jj2lsj/source-manifest.json).
Source and numerical behavior are described in [the implementation guide](../../docs/implemented/jj2lsj_openmp.md).

From the suite root:

```sh
cmake --build build --target jj2lsj_omp test.jj2lsj_coefficients --parallel 4
ctest --test-dir build -R '^jj2lsj_' --output-on-failure
```

Use `-DGRASP_JJ2LSJ_OPENMP=OFF` in a separate build to test the serial-only
configuration. Debug validation uses `-DCMAKE_Fortran_FLAGS_DEBUG='-g -fcheck=all -fbacktrace'`.
These tests establish numerical and file-format compatibility for the fixtures;
they are not a production performance benchmark or an atomic-accuracy study.


The postprocessing regression also checks optional `JJ2LSJ_PROFILE=1`, exact
reported output sizes, serial versus parallel ASF ranking, and unchanged unique
labels. `test.jj2lsj_ranking` compares against an independent insertion-sort
oracle, including signed zeros, equal squares, nonexact sums and an 84,000-row
zero plateau. Nonfinite coefficients must fail explicitly. These tests preserve
the original LS-order norm accumulation and all original eigenvector entries.

`JJ2LSJ_POST_SERIAL=1` disables only ranking concurrency; `JJ2LSJ_SERIAL=1`
disables both conversion and ranking concurrency. Both use the new stable sort.
The original executable is needed to compare against the old sorting algorithm.
