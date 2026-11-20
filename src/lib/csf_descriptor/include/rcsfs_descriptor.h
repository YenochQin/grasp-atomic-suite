#ifndef RCSFS_DESCRIPTOR_H
#define RCSFS_DESCRIPTOR_H
#include <stdint.h>
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
/* ABI v1: status 0 is success. Errors write a NUL-terminated message if
 * capacity > 0. Handles are exclusive, independent, and closed by address.
 * Initialize an open() handle slot to NULL; close before reopening. A live
 * slot is rejected unchanged. Paths are UTF-8, NUL terminated.
 * No caller buffer pointer is retained.
 * Info: core, peel, blocks, rows, electrons (five int64_t entries).
 * Tables: orbital (n,kappa) pairs, block lengths, block (2J,parity) pairs.
 * Table capacities count orbital/block records, not individual scalar entries.
 * Read buffer: channel + row_width*row. On failure output is invalid.
 * All pointers and sizes must describe valid caller-owned storage. */
int32_t rcsfs_descriptor_abi_version(void);
/* Header-only radial input. Orbitals are (n,kappa) pairs in core/peel order.
 * Info contains orbital count, core count, and total electrons. */
int32_t rcsfs_descriptor_read_radial_header(const char*, int32_t*, size_t, int32_t*, char*, size_t);
int32_t rcsfs_descriptor_open(const char*, const char*, void**, char*, size_t);
int32_t rcsfs_descriptor_info(void*, int64_t*, char*, size_t);
int32_t rcsfs_descriptor_tables(void*, int32_t*, size_t, int64_t*, int32_t*, size_t, char*, size_t);
/* One count per block. Returns an error if the header has no partition. */
int32_t rcsfs_descriptor_zero_order_counts(void*, int64_t*, size_t, char*, size_t);
int32_t rcsfs_descriptor_read(void*, int8_t*, size_t, size_t, char*, size_t);
void rcsfs_descriptor_close(void**);
#ifdef __cplusplus
}
#endif
#endif
