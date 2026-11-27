//! C ABI v3. Callers must provide valid, aligned buffers for the stated sizes.
//! Output buffers are caller-owned, exclusively writable during the call, and
//! disjoint from handle-owned storage.
//! An error buffer must be writable for its capacity when non-null and capacity > 0.
//! Calls using the same live handle must be serialized; separate handles are independent.
//! The unmangled `rcsfs_descriptor_*` symbols are reserved for this ABI and must
//! have a single definition in the linked program.
//!
//! Raw pointers are necessary for the Fortran C ABI. The descriptor reader
//! itself takes Rust references and returns ordinary Rust values and errors.
use crate::descriptor_reader::{DescriptorReader, RadialHeader};
use anyhow::{Result, ensure};
use std::ffi::{CStr, c_char, c_void};
use std::panic::{self, AssertUnwindSafe};
use std::path::Path;
use std::ptr;

struct Handle {
    reader: DescriptorReader,
    batch: Vec<i8>,
    offset: usize,
}

/// # Safety
/// When non-null and `capacity > 0`, `error` must remain writable for `capacity`
/// bytes after `action` finishes.
unsafe fn invoke(error: *mut c_char, capacity: usize, action: impl FnOnce() -> Result<()>) -> i32 {
    // Recoverable input and I/O failures use Result. This boundary also
    // reports an unexpected unwind before it can reach the Fortran caller.
    let result = panic::catch_unwind(AssertUnwindSafe(action));
    let message = match result {
        Ok(Ok(())) => String::new(),
        Ok(Err(e)) => format!("{e:#}"),
        Err(_) => "Rust descriptor reader panicked".into(),
    };
    if !error.is_null() && capacity > 0 {
        let bytes = message.as_bytes();
        let mut length = bytes.len().min(capacity - 1);
        while !message.is_char_boundary(length) {
            length -= 1;
        }
        // SAFETY: The caller provides writable storage disjoint from the local
        // message. The copy and trailing NUL fit because length < capacity.
        unsafe {
            ptr::copy_nonoverlapping(bytes.as_ptr(), error.cast(), length);
            *error.add(length) = 0;
        }
    }
    if message.is_empty() { 0 } else { 1 }
}

/// The C ABI revision implemented by the exported descriptor functions.
#[unsafe(no_mangle)]
pub extern "C" fn rcsfs_descriptor_abi_version() -> i32 {
    3
}

/// The on-disk format accepted and returned by this reader.
#[unsafe(no_mangle)]
pub extern "C" fn rcsfs_descriptor_format_version() -> i32 {
    i32::from(crate::descriptor_schema::DESCRIPTOR_VERSION)
}

/// Copy the verified Parquet SHA-256 as 64 ASCII bytes. Its metadata binds the ordered header.
/// # Safety
/// The handle must be live and the output must be writable for 64 bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rcsfs_descriptor_basis_id(handle: *const c_void, output: *mut c_char) {
    // SAFETY: The caller provides a live handle and disjoint writable storage.
    unsafe {
        let reader = &(*handle.cast::<Handle>()).reader;
        ptr::copy_nonoverlapping(
            reader.basis_id().as_ptr(),
            output.cast(),
            reader.basis_id().len(),
        );
    }
}

/// Read radial input metadata from TOML. No descriptor file or handle is used.
///
/// # Safety
/// The path must be readable through its NUL terminator and stay unchanged while
/// borrowed. The caller must supply writable storage for `2 * orbital_capacity`
/// i32 orbital entries and three i32 info entries.
/// The error buffer, if non-null, must be writable for `capacity` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rcsfs_descriptor_read_radial_header(
    header: *const c_char,
    orbitals: *mut i32,
    orbital_capacity: usize,
    info: *mut i32,
    error: *mut c_char,
    capacity: usize,
) -> i32 {
    let action = || {
        ensure!(
            !header.is_null() && !orbitals.is_null() && !info.is_null(),
            "null radial header argument"
        );
        // SAFETY: The checked pointer refers to a caller-provided C string
        // readable through its NUL terminator and unchanged while borrowed.
        let header_path = unsafe { CStr::from_ptr(header) };
        let header = RadialHeader::read(Path::new(header_path.to_str()?))?;
        ensure!(
            header.orbitals.len() <= orbital_capacity,
            "orbital table exceeds configured max_orbitals capacity"
        );
        let values = [
            i32::try_from(header.orbitals.len())?,
            i32::try_from(header.core_count)?,
            header.electron_count,
        ];
        for (i, orbital) in header.orbitals.into_iter().enumerate() {
            // SAFETY: The caller provides room for two i32s per orbital, and
            // the capacity check bounds both indices.
            unsafe {
                *orbitals.add(2 * i) = orbital.principal_quantum_number;
                *orbitals.add(2 * i + 1) = orbital.kappa;
            }
        }
        // SAFETY: info points to three writable i32s disjoint from local values.
        unsafe { ptr::copy_nonoverlapping(values.as_ptr(), info, values.len()) };
        Ok(())
    };
    // SAFETY: The caller provides an optional error buffer valid for capacity bytes.
    unsafe { invoke(error, capacity, action) }
}

/// # Safety
/// Both paths must be readable through their NUL terminators and stay unchanged
/// while borrowed. Output must point to an initialized, writable handle slot
/// with exclusive access for this call. Initialize it to NULL before opening.
/// A live slot is rejected unchanged; close it before opening another file.
/// The optional error buffer follows the module's capacity contract.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rcsfs_descriptor_open(
    path: *const c_char,
    header: *const c_char,
    output: *mut *mut c_void,
    error: *mut c_char,
    capacity: usize,
) -> i32 {
    let action = || {
        ensure!(!output.is_null(), "null output handle");
        // SAFETY: The caller provides a valid, initialized handle slot.
        let current = unsafe { *output };
        ensure!(
            current.is_null(),
            "close the existing handle before reopening"
        );
        ensure!(!path.is_null() && !header.is_null(), "null input path");
        // SAFETY: The caller's path is readable through its NUL terminator
        // and unchanged while borrowed.
        let path = unsafe { CStr::from_ptr(path) }.to_str()?;
        // SAFETY: The caller's header path has the same C string guarantees.
        let header = unsafe { CStr::from_ptr(header) }.to_str()?;
        let reader = DescriptorReader::open(Path::new(path), Path::new(header), 8192)?;
        let handle = Box::into_raw(Box::new(Handle {
            reader,
            batch: Vec::new(),
            offset: 0,
        }));
        // SAFETY: The caller grants write access to the slot, which was NULL.
        unsafe { *output = handle.cast() };
        Ok(())
    };
    // SAFETY: The caller provides an optional error buffer valid for capacity bytes.
    unsafe { invoke(error, capacity, action) }
}

/// # Safety
/// Handle must be live and used serially. The optional error buffer follows
/// the module's capacity contract.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rcsfs_descriptor_check_spin_angular(
    handle: *mut c_void,
    error: *mut c_char,
    capacity: usize,
) -> i32 {
    let action = || {
        ensure!(!handle.is_null(), "null descriptor handle");
        // SAFETY: the caller supplies a live, exclusively accessed handle from open.
        let h = unsafe { &mut *handle.cast::<Handle>() };
        h.reader.check_spin_angular()
    };
    // SAFETY: the caller supplies a writable error buffer for capacity bytes.
    unsafe { invoke(error, capacity, action) }
}

/// # Safety
/// Handle must be live and used serially; info has five writable i64 entries.
/// The optional error buffer follows the module's capacity contract.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rcsfs_descriptor_info(
    handle: *mut c_void,
    info: *mut i64,
    error: *mut c_char,
    capacity: usize,
) -> i32 {
    let action = || {
        ensure!(!handle.is_null() && !info.is_null(), "null info argument");
        // SAFETY: The handle came from open and stays live with no concurrent
        // access for the duration of this call.
        let reader = unsafe { &(*handle.cast::<Handle>()).reader };
        let values = [
            i64::try_from(reader.core_count())?,
            i64::try_from(reader.layout().subshell_count)?,
            i64::try_from(reader.block_lengths().len())?,
            i64::try_from(reader.block_lengths().iter().sum::<u64>())?,
            i64::from(reader.electron_count()),
        ];
        // SAFETY: info points to five writable i64s disjoint from local values.
        unsafe { ptr::copy_nonoverlapping(values.as_ptr(), info, 5) };
        Ok(())
    };
    // SAFETY: The caller provides an optional error buffer valid for capacity bytes.
    unsafe { invoke(error, capacity, action) }
}

/// # Safety
/// Handle must be live and used serially. The caller provides writable arrays
/// for the stated capacities. Orbital and symmetry arrays are interleaved pairs,
/// and lengths has one entry per block. The optional error buffer follows the
/// module's capacity contract.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rcsfs_descriptor_tables(
    handle: *mut c_void,
    orbitals: *mut i32,
    orbital_capacity: usize,
    lengths: *mut i64,
    symmetries: *mut i32,
    block_capacity: usize,
    error: *mut c_char,
    capacity: usize,
) -> i32 {
    let action = || {
        ensure!(
            !handle.is_null() && !orbitals.is_null() && !lengths.is_null() && !symmetries.is_null(),
            "null tables argument"
        );
        // SAFETY: The handle came from open and stays live with no concurrent
        // access for the duration of this call.
        let reader = unsafe { &(*handle.cast::<Handle>()).reader };
        ensure!(
            orbital_capacity >= reader.orbitals().len()
                && block_capacity >= reader.block_lengths().len(),
            "table buffer too small"
        );
        for (i, orbital) in reader.orbitals().iter().enumerate() {
            // SAFETY: The caller provides two writable i32s per orbital, and
            // the capacity check bounds both indices.
            unsafe {
                *orbitals.add(2 * i) = orbital.principal_quantum_number;
                *orbitals.add(2 * i + 1) = orbital.kappa;
            }
        }
        let blocks = reader.block_lengths().iter().zip(reader.block_symmetries());
        for (i, (&length, symmetry)) in blocks.enumerate() {
            let length = i64::try_from(length)?;
            // SAFETY: The caller provides one length and two symmetry entries
            // per block. The loop and capacity check bound all three indices.
            unsafe {
                *lengths.add(i) = length;
                *symmetries.add(2 * i) = symmetry.total_two_j;
                *symmetries.add(2 * i + 1) = symmetry.parity;
            }
        }
        Ok(())
    };
    // SAFETY: The caller provides an optional error buffer valid for capacity bytes.
    unsafe { invoke(error, capacity, action) }
}

/// # Safety
/// Handle must be live and used serially. Counts has `block_capacity` writable
/// i64 entries. The optional error buffer follows the module's capacity contract.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rcsfs_descriptor_zero_order_counts(
    handle: *mut c_void,
    counts: *mut i64,
    block_capacity: usize,
    error: *mut c_char,
    capacity: usize,
) -> i32 {
    let action = || {
        ensure!(
            !handle.is_null() && !counts.is_null(),
            "null zero-order counts argument"
        );
        // SAFETY: The handle came from open and stays live with no concurrent
        // access for the duration of this call.
        let reader = unsafe { &(*handle.cast::<Handle>()).reader };
        let values = reader.zero_order_counts().ok_or_else(|| {
            anyhow::anyhow!(
                "missing block_info.zero_order_counts; regenerate the dataset with rcsfs zero-first"
            )
        })?;
        ensure!(
            block_capacity >= values.len(),
            "zero-order counts buffer too small"
        );
        // SAFETY: The capacity check ensures the copy fits in the caller's
        // writable counts buffer, which is disjoint from the reader's values.
        unsafe { ptr::copy_nonoverlapping(values.as_ptr(), counts, values.len()) };
        Ok(())
    };
    // SAFETY: The caller provides an optional error buffer valid for capacity bytes.
    unsafe { invoke(error, capacity, action) }
}

/// Select a zero-based explicit block. Discard any buffered rows from an earlier selection.
///
/// # Safety
/// Handle must be live and exclusively borrowed. The optional error buffer follows
/// the module's capacity contract.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rcsfs_descriptor_select_block(
    handle: *mut c_void,
    block: usize,
    error: *mut c_char,
    capacity: usize,
) -> i32 {
    let action = || {
        ensure!(!handle.is_null(), "null descriptor handle");
        // SAFETY: the caller provides exclusive access to a live handle.
        let h = unsafe { &mut *handle.cast::<Handle>() };
        h.reader.select_block(block)?;
        h.batch.clear();
        h.offset = 0;
        Ok(())
    };
    // SAFETY: the caller supplies a writable error buffer for capacity bytes.
    unsafe { invoke(error, capacity, action) }
}

/// # Safety
/// Handle must be live and exclusively borrowed. Output has `size` writable
/// bytes. The optional error buffer follows the module's capacity contract.
/// Reads exactly rows; a failure invalidates the output buffer.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rcsfs_descriptor_read(
    handle: *mut c_void,
    output: *mut i8,
    size: usize,
    rows: usize,
    error: *mut c_char,
    capacity: usize,
) -> i32 {
    let action = || {
        ensure!(!handle.is_null() && !output.is_null(), "null read argument");
        // SAFETY: The handle came from open and the caller grants exclusive
        // access to the live Handle for this call.
        let h = unsafe { &mut *handle.cast::<Handle>() };
        let count = rows
            .checked_mul(h.reader.layout().descriptor_row_channel_count())
            .ok_or_else(|| anyhow::anyhow!("size overflow"))?;
        ensure!(size >= count, "output buffer too small");
        let mut written = 0;
        while written < count {
            if h.offset == h.batch.len() {
                // Keep the handle's batch empty until the reader validates the
                // replacement. Reuse the allocation on success.
                let mut batch = std::mem::take(&mut h.batch);
                h.offset = 0;
                ensure!(h.reader.read_batch(&mut batch)? > 0, "unexpected EOF");
                h.batch = batch;
            }
            let source = &h.batch[h.offset..];
            let n = (count - written).min(source.len());
            // SAFETY: n fits in the source slice and in the remaining checked
            // output span. The caller's output is disjoint from the batch.
            unsafe { ptr::copy_nonoverlapping(source.as_ptr(), output.add(written), n) };
            written += n;
            h.offset += n;
        }
        Ok(())
    };
    // SAFETY: The caller provides an optional error buffer valid for capacity bytes.
    unsafe { invoke(error, capacity, action) }
}

/// # Safety
/// A non-null pointer must address an initialized, writable slot containing
/// NULL or a live handle returned by open. The caller must have exclusive
/// access to the slot and handle, and must not use any copies after closing.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rcsfs_descriptor_close(handle: *mut *mut c_void) {
    if !handle.is_null() {
        // SAFETY: The caller provides an initialized handle slot.
        let current = unsafe { *handle };
        if !current.is_null() {
            // SAFETY: The caller hands back ownership of the live Box from
            // open and grants write access to the slot so it can be cleared.
            unsafe {
                drop(Box::from_raw(current.cast::<Handle>()));
                *handle = ptr::null_mut();
            }
        }
    }
}
