use super::*;
use arrow_array::RecordBatch;
use parquet::{arrow::ArrowWriter, file::properties::WriterProperties};
use std::sync::Arc;

fn fixture(
    root: &Path,
    mutate: impl FnOnce(&mut Vec<parquet::file::metadata::KeyValue>),
) -> Result<(std::path::PathBuf, std::path::PathBuf)> {
    let header = root.join("basis_header.toml");
    let path = root.join("basis_descriptors.parquet");
    std::fs::write(
        &header,
        "[header_info]\nheader_lines = ['Core subshells:', '1s', 'Peel subshells:', '2s 2p', 'CSF(s):']\nelectron_count=4\n[block_info]\nblock_count=2\nblock_lengths=[1,1]\n",
    )?;
    let layout = DescriptorLayout::new(DescriptorVersion::V2, 2);
    let schema = crate::descriptor_schema::output_schema(layout)?;
    let hash = crate::descriptor_schema::hash_header_file(&header)?;
    let mut metadata = crate::descriptor_schema::output_kv_metadata(
        layout,
        &["2s".into(), "2p".into()],
        Some(&hash),
        None,
    );
    mutate(&mut metadata);
    let properties = WriterProperties::builder()
        .set_key_value_metadata(Some(metadata))
        .set_max_row_group_row_count(Some(1))
        .build();
    let mut writer = ArrowWriter::try_new(File::create(&path)?, schema.clone(), Some(properties))?;
    let rows = [
        [2, -1, 0, -1, 0, -1, 0, -1, 0, 1],
        [0, -1, 0, -1, 2, 0, 0, -1, 0, 1],
    ];
    let columns = (0..10)
        .map(|c| Arc::new(Int8Array::from(vec![rows[0][c], rows[1][c]])) as arrow_array::ArrayRef)
        .collect();
    writer.write(&RecordBatch::try_new(schema, columns)?)?;
    writer.close()?;
    Ok((path, header))
}

#[test]
fn strict_footer_keys_and_ordered_batches() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let required = [
        "descriptor_version",
        "descriptor_dtype",
        "max_two_j",
        "seniority_encoding",
        "channels_per_subshell",
        "subshell_count",
        "missing_sentinel",
        "normalized",
        "source_header_sha256",
        "peel_subshells",
        "feature_columns",
        "global_columns",
    ];
    for key in required {
        let (path, header) = fixture(directory.path(), |metadata| {
            metadata.retain(|entry| entry.key != key)
        })?;
        assert!(
            DescriptorReader::open(&path, &header, 1).is_err(),
            "missing {key}"
        );
        let (path, header) = fixture(directory.path(), |metadata| {
            let duplicate = metadata
                .iter()
                .find(|entry| entry.key == key)
                .unwrap()
                .clone();
            metadata.push(duplicate);
        })?;
        assert!(
            DescriptorReader::open(&path, &header, 1).is_err(),
            "duplicate {key}"
        );
        let (path, header) = fixture(directory.path(), |metadata| {
            metadata
                .iter_mut()
                .find(|entry| entry.key == key)
                .unwrap()
                .value = Some("invalid".into());
        })?;
        assert!(
            DescriptorReader::open(&path, &header, 1).is_err(),
            "invalid {key}"
        );
    }
    let (path, header) = fixture(directory.path(), |_| {})?;
    let mut first = DescriptorReader::open(&path, &header, 1)?;
    let mut second = DescriptorReader::open(&path, &header, 2)?;
    assert_eq!(first.block_lengths, [1, 1]);
    assert_eq!(first.block_symmetries, [(0, 1), (0, 1)]);
    let (mut a, mut b) = (Vec::new(), Vec::new());
    assert_eq!(first.read_batch(&mut a)?, 1);
    assert_eq!(second.read_batch(&mut b)?, 2);
    assert_eq!(a, b[..10]);
    assert_eq!(first.read_batch(&mut a)?, 1);
    assert_eq!(a, b[10..]);
    assert_eq!(first.read_batch(&mut a)?, 0);
    Ok(())
}

#[test]
fn rejects_corrupt_payload_and_header_semantics() -> Result<()> {
    let directory = tempfile::tempdir()?;
    for replacement in [
        "block_lengths=[0,2]",
        "block_lengths=[1]",
        "block_lengths=[-1,3]",
    ] {
        let (path, header) = fixture(directory.path(), |_| {})?;
        let text = std::fs::read_to_string(&header)?.replace("block_lengths=[1,1]", replacement);
        std::fs::write(&header, text)?;
        assert!(DescriptorReader::open(&path, &header, 1).is_err());
    }
    let (path, header) = fixture(directory.path(), |_| {})?;
    let original = std::fs::read(&path)?;
    for bytes in [
        &original[..original.len() - 1],
        &original[..original.len() / 2],
        b"GCSDESC1",
    ] {
        std::fs::write(&path, bytes)?;
        assert!(DescriptorReader::open(&path, &header, 1).is_err());
    }
    Ok(())
}

#[cfg(feature = "descriptor-ffi")]
#[test]
fn ffi_buffers_errors_and_independent_handles() -> Result<()> {
    use crate::descriptor_ffi::*;
    use std::ffi::{CString, c_char, c_void};
    let directory = tempfile::tempdir()?;
    let (path, header) = fixture(directory.path(), |_| {})?;
    let path = CString::new(path.to_str().unwrap())?;
    let header = CString::new(header.to_str().unwrap())?;
    let mut handles = [std::ptr::null_mut::<c_void>(); 2];
    let mut error = [0 as c_char; 512];
    unsafe {
        assert_eq!(rcsfs_descriptor_abi_version(), 1);
        for handle in &mut handles {
            assert_eq!(
                rcsfs_descriptor_open(
                    path.as_ptr(),
                    header.as_ptr(),
                    handle,
                    error.as_mut_ptr(),
                    512
                ),
                0
            );
        }
        let original = handles[0];
        assert_ne!(
            rcsfs_descriptor_open(
                path.as_ptr(),
                header.as_ptr(),
                &mut handles[0],
                error.as_mut_ptr(),
                512
            ),
            0
        );
        assert_eq!(handles[0], original);
        let mut info = [0i64; 5];
        assert_eq!(
            rcsfs_descriptor_info(handles[0], info.as_mut_ptr(), error.as_mut_ptr(), 512),
            0
        );
        assert_eq!(info, [1, 2, 2, 2, 4]);
        let mut orbitals = [-99i32; 6];
        let mut lengths = [0i64; 2];
        let mut symmetries = [0i32; 4];
        assert_ne!(
            rcsfs_descriptor_tables(
                handles[0],
                orbitals.as_mut_ptr(),
                1,
                lengths.as_mut_ptr(),
                symmetries.as_mut_ptr(),
                2,
                error.as_mut_ptr(),
                512
            ),
            0
        );
        assert_eq!(orbitals, [-99; 6]);
        assert_eq!(
            rcsfs_descriptor_tables(
                handles[0],
                orbitals.as_mut_ptr(),
                3,
                lengths.as_mut_ptr(),
                symmetries.as_mut_ptr(),
                2,
                error.as_mut_ptr(),
                512
            ),
            0
        );
        assert_eq!(orbitals, [1, -1, 2, -1, 2, -2]);
        assert_eq!(lengths, [1, 1]);
        assert_eq!(symmetries, [0, 1, 0, 1]);
        let mut rows = [-99i8; 20];
        assert_ne!(
            rcsfs_descriptor_read(handles[0], rows.as_mut_ptr(), 1, 1, error.as_mut_ptr(), 512),
            0
        );
        assert_eq!(rows, [-99; 20]);
        assert_eq!(
            rcsfs_descriptor_read(
                handles[0],
                rows.as_mut_ptr(),
                20,
                2,
                error.as_mut_ptr(),
                512
            ),
            0
        );
        assert_eq!(
            rows,
            [
                2, -1, 0, -1, 0, -1, 0, -1, 0, 1, 0, -1, 0, -1, 2, 0, 0, -1, 0, 1
            ]
        );
        let mut other = [0i8; 10];
        assert_eq!(
            rcsfs_descriptor_read(
                handles[1],
                other.as_mut_ptr(),
                10,
                1,
                error.as_mut_ptr(),
                512
            ),
            0
        );
        assert_eq!(other, rows[..10]);
        assert_ne!(
            rcsfs_descriptor_read(
                handles[0],
                other.as_mut_ptr(),
                10,
                1,
                error.as_mut_ptr(),
                512
            ),
            0
        );
        for handle in &mut handles {
            rcsfs_descriptor_close(handle);
            rcsfs_descriptor_close(handle);
            assert!(handle.is_null());
        }
        assert_ne!(
            rcsfs_descriptor_info(handles[0], info.as_mut_ptr(), error.as_mut_ptr(), 512),
            0
        );
    }
    Ok(())
}

#[test]
fn closed_core_batches_and_header_binding() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let header = dir.path().join("header.toml");
    let path = dir.path().join("descriptors.parquet");
    let text = "[header_info]\nheader_lines = ['Core subshells:', '1s', 'Peel subshells:', '', 'CSF(s):']\n[block_info]\nblock_lengths = [1]\n";
    std::fs::write(&header, text)?;
    let layout = DescriptorLayout::new(DescriptorVersion::V2, 0);
    let schema = crate::descriptor_schema::output_schema(layout)?;
    let hash = crate::descriptor_schema::hash_header_file(&header)?;
    let kv = crate::descriptor_schema::output_kv_metadata(layout, &[], Some(&hash), None);
    let properties = WriterProperties::builder()
        .set_key_value_metadata(Some(kv))
        .build();
    let mut writer = ArrowWriter::try_new(File::create(&path)?, schema.clone(), Some(properties))?;
    writer.write(&RecordBatch::try_new(
        schema,
        vec![
            Arc::new(Int8Array::from(vec![0])),
            Arc::new(Int8Array::from(vec![1])),
        ],
    )?)?;
    writer.close()?;
    let mut reader = DescriptorReader::open(&path, &header, 1)?;
    let mut values = Vec::new();
    assert_eq!(reader.read_batch(&mut values)?, 1);
    assert_eq!(values, [0, 1]);
    assert_eq!(reader.read_batch(&mut values)?, 0);
    assert!(values.is_empty());
    std::fs::write(&header, format!("{text}\n"))?;
    assert!(DescriptorReader::open(&path, &header, 1).is_err());
    assert!(DescriptorReader::open(&path, &header, 0).is_err());
    Ok(())
}
