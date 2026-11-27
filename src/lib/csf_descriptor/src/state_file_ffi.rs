//! Private file-integrity and locking operations for the Suite Fortran applications.
//!
//! The Fortran C ABI passes paths and outputs as raw pointers and transfers file
//! ownership as a POSIX descriptor. Unsafe operations are needed to translate
//! those ABI values to Rust borrows and an owned file. Rust's safe I/O and
//! locking APIs perform the file operations.

use sha2::{Digest, Sha256};
use std::ffi::{CStr, OsStr, c_char};
use std::fs::{File, OpenOptions};
use std::io;
use std::os::fd::{FromRawFd, IntoRawFd, RawFd};
use std::os::unix::ffi::OsStrExt;

/// Write the file's SHA-256 as 64 lowercase ASCII bytes and return zero on success.
///
/// # Safety
///
/// `path` must be readable through its NUL terminator and remain unchanged
/// during the call. `digest` must point to 64 writable bytes, disjoint from
/// `path` and exclusively accessible for the write.
#[unsafe(no_mangle)]
unsafe extern "C" fn grasp_file_sha256(path: *const c_char, digest: *mut u8) -> i32 {
    // SAFETY: The caller provides a readable, immutable C string for this call.
    let path = OsStr::from_bytes(unsafe { CStr::from_ptr(path) }.to_bytes());
    let Ok(mut file) = File::open(path) else {
        return 1;
    };
    let mut hash = Sha256::new();
    if io::copy(&mut file, &mut hash).is_err() {
        return 1;
    }
    let encoded = format!("{:x}", hash.finalize());
    // SAFETY: SHA-256 formats to 64 hex bytes in local storage. The caller
    // provides 64 writable bytes disjoint from that storage.
    unsafe { digest.copy_from_nonoverlapping(encoded.as_ptr(), 64) };
    0
}

/// Acquire a shared lock for zero or an exclusive lock for a nonzero `writing`.
/// Return an owned descriptor on success or `-1` on failure.
///
/// # Safety
///
/// `path` must be readable through its NUL terminator and remain unchanged
/// during the call. The caller takes ownership of a successful result and
/// must return it once to [`grasp_unlock_state_files`].
#[unsafe(no_mangle)]
unsafe extern "C" fn grasp_lock_state_files(path: *const c_char, writing: i32) -> RawFd {
    // SAFETY: The caller provides a readable, immutable C string for this call.
    let path = OsStr::from_bytes(unsafe { CStr::from_ptr(path) }.to_bytes());
    let writing = writing != 0;
    let Ok(file) = OpenOptions::new()
        .read(true)
        .write(writing)
        .create(writing)
        .truncate(false)
        .open(path)
    else {
        return -1;
    };
    let lock_result = if writing {
        file.lock()
    } else {
        file.lock_shared()
    };
    if lock_result.is_err() {
        return -1;
    }
    file.into_raw_fd()
}

/// Close the lock file and release its lock.
///
/// # Safety
///
/// `lock` must be the still-open descriptor returned by a successful
/// [`grasp_lock_state_files`] call. The caller relinquishes ownership and must
/// neither use nor close it after this call.
#[unsafe(no_mangle)]
unsafe extern "C" fn grasp_unlock_state_files(lock: RawFd) {
    // SAFETY: The caller transfers exclusive ownership of a live file descriptor.
    drop(unsafe { File::from_raw_fd(lock) });
}
