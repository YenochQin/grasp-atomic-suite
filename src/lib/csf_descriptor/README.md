# Suite CSF descriptor reader

This directory owns the suite's native Parquet/header reader. CMake builds it
from these local sources and places `libcsf_descriptor.a` beside the other
suite libraries in `<build>/lib/`; installation copies it to `<prefix>/lib/`.
Neither an rCSFs checkout nor a prebuilt rCSFs library is used by this build.

The initial implementation was copied from
[rCSFs](https://github.com/YenochQin/rCSFs), revision
`f17a862d77c45164a6d31ae1e9cfeb4efbd04cf6` (2026-11-17):

- The descriptor transport, reader, and C interface originate from that
  revision. The suite reader's tests are in `src/tests.rs`.
- `include/rcsfs_descriptor.h` declares ABI 3 and exposes the accepted format
  through `rcsfs_descriptor_format_version`. The separate radial header entry point serves `rwfnestimate`.
  It shares the TOML parser and returns orbital metadata from the header alone.
- The manifest and lockfile originate from `native-reader/`; only the local
  package/library identity, source entry point, and library kinds are changed.
- The original MIT copyright and permission notice are retained in `LICENSE`.

Future reader fixes belong here and are tested with the suite. Upstream changes
must be reviewed before being copied; building never synchronizes external
source. The reader accepts indexed V3. Opening verifies the full Parquet hash and bound
block directory. Block selection restricts decoding to the indexed rows. Rust
validates the subshell states and cumulative couplings as rows are read. The Fortran adapter installs the validated rows in the
Suite's orbital order.

Rust 1.99.0 or a newer release, Cargo, and a C compiler are needed at build time.
The lockfile pins third-party dependencies. CMake runs Cargo with `--locked`
and allows Cargo to download missing dependencies. Set
`-DGRASP_CSF_READER_OFFLINE=ON` to use the local dependency cache in offline mode.
This does not fetch rCSFs. Runtime needs neither Python nor Rust tools.

The Rust format, FFI, and lifecycle tests are registered as
`csf_descriptor_rust` in CTest. They can also be run independently, from the
suite root:

```fish
cargo test --manifest-path src/lib/csf_descriptor/Cargo.toml \
    --target-dir build-reader-tests --release --locked --offline --jobs 4
```

The default and extended reader contracts are documented in
[native CSF inputs](../../../docs/implemented/native_csf_descriptors.md).
