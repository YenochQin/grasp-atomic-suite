#ifndef RCSFS_DESCRIPTOR_H
#define RCSFS_DESCRIPTOR_H
#include <stdint.h>
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
/* ABI v2: status 0 is success. Errors write a NUL-terminated message if
 * capacity > 0. Handles are exclusive, independent, and closed by address.
 * Initialize an open() handle slot to NULL; close before reopening. A live
 * slot is rejected unchanged. Paths are UTF-8, NUL terminated.
 * No caller buffer pointer is retained.
 * Info: core, peel, blocks, rows, electrons (five int64_t entries).
 * Tables: orbital (n,kappa) pairs, block lengths, block (2J,parity) pairs.
 * Table capacities count orbital/block records, not individual scalar entries.
 * Read buffer: V3 codes at channel + row_width*row. On failure output is invalid.
 * All pointers and sizes must describe valid caller-owned storage. */
int32_t rcsfs_descriptor_abi_version(void);
int32_t rcsfs_descriptor_format_version(void);
/* Header-only radial input. Orbitals are (n,kappa) pairs in core/peel order.
 * Info contains orbital count, core count, and total electrons. */
int32_t rcsfs_descriptor_read_radial_header(
    const char *header_path,
    int32_t *orbital_pairs,
    size_t orbital_capacity,
    int32_t *radial_header_info,
    char *error_buffer,
    size_t error_buffer_capacity);
int32_t rcsfs_descriptor_open(
    const char *descriptor_path,
    const char *header_path,
    void **output_handle,
    char *error_buffer,
    size_t error_buffer_capacity);
/* Check the spin-angular kernel's occupation limits after the opening scan. */
int32_t rcsfs_descriptor_check_spin_angular(void *handle, char *error_buffer, size_t error_buffer_capacity);
int32_t rcsfs_descriptor_info(
    void *handle,
    int64_t *descriptor_info,
    char *error_buffer,
    size_t error_buffer_capacity);
int32_t rcsfs_descriptor_tables(
    void *handle,
    int32_t *orbital_pairs,
    size_t orbital_capacity,
    int64_t *block_lengths,
    int32_t *block_symmetry_pairs,
    size_t block_capacity,
    char *error_buffer,
    size_t error_buffer_capacity);
/* One count per block. Returns an error if the header has no partition. */
int32_t rcsfs_descriptor_zero_order_counts(
    void *handle,
    int64_t *zero_order_counts,
    size_t block_capacity,
    char *error_buffer,
    size_t error_buffer_capacity);
int32_t rcsfs_descriptor_read(
    void *handle,
    int8_t *descriptor_rows,
    size_t output_capacity,
    size_t row_count,
    char *error_buffer,
    size_t error_buffer_capacity);
void rcsfs_descriptor_close(void **handle);
#ifdef __cplusplus
}
#endif
#endif
