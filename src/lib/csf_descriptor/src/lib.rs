//! Suite-owned native CSF descriptor reader; no Python or generator dependency.
#[cfg(feature = "descriptor-ffi")]
pub mod descriptor_ffi;
pub mod descriptor_reader;
#[path = "descriptor_transport.rs"]
pub mod descriptor_schema;
