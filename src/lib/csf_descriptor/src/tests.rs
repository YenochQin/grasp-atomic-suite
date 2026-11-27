use super::{BlockSymmetry, DescriptorReader};
use crate::descriptor_schema::DescriptorLayout;
use anyhow::Result;
use arrow_array::{Int8Array, RecordBatch};
use parquet::arrow::ArrowWriter;
use parquet::file::metadata::KeyValue;
use parquet::file::properties::WriterProperties;
use sha2::{Digest, Sha256};
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

fn fixture_metadata(
    layout: DescriptorLayout,
    peel: &[&str],
    header: &Path,
) -> Result<Vec<KeyValue>> {
    Ok([
        ("descriptor_version", "3".into()),
        ("descriptor_dtype", "int8".into()),
        ("max_two_j", "126".into()),
        ("seniority_encoding", "physical".into()),
        ("channels_per_subshell", "4".into()),
        ("subshell_count", peel.len().to_string()),
        ("normalized", "false".into()),
        (
            "source_header_sha256",
            format!("{:x}", Sha256::digest(std::fs::read(header)?)),
        ),
        ("peel_subshells", serde_json::to_string(peel)?),
        (
            "feature_columns",
            serde_json::to_string(&layout.feature_column_names())?,
        ),
        (
            "global_columns",
            serde_json::to_string(&layout.global_column_names())?,
        ),
    ]
    .into_iter()
    .map(|(key, value)| KeyValue::new(key.into(), Some(value)))
    .collect())
}

fn seal_fixture(path: &Path, header: &Path) -> Result<()> {
    let hash = format!("{:x}", Sha256::digest(std::fs::read(path)?));
    let body = std::fs::read_to_string(header)?;
    std::fs::write(header, format!("descriptor_sha256 = \"{hash}\"\n{body}"))?;
    Ok(())
}

fn fixture(
    root: &Path,
    mutate: impl FnOnce(&mut Vec<parquet::file::metadata::KeyValue>),
) -> Result<(std::path::PathBuf, std::path::PathBuf)> {
    fixture_with_index(
        root,
        "[{two_j=0,parity=1,start_row=0,length=1},{two_j=0,parity=1,start_row=1,length=1}]",
        1,
        mutate,
    )
}

fn fixture_with_index(
    root: &Path,
    blocks: &str,
    group_size: usize,
    mutate: impl FnOnce(&mut Vec<KeyValue>),
) -> Result<(std::path::PathBuf, std::path::PathBuf)> {
    let header = root.join("basis_header.toml");
    let path = root.join("basis_descriptors.parquet");
    std::fs::write(
        &header,
        format!(
            "[header_info]\nheader_lines = ['Core subshells:', '1s', 'Peel subshells:', '2s 2p', 'CSF(s):']\nelectron_count=4\n[block_info]\nblock_count=2\nblock_lengths=[1,1]\nblocks={blocks}\n"
        ),
    )?;
    let layout = DescriptorLayout { subshell_count: 2 };
    let schema = Arc::new(layout.schema());
    let mut metadata = fixture_metadata(layout, &["2s", "2p"], &header)?;
    mutate(&mut metadata);
    let properties = WriterProperties::builder()
        .set_key_value_metadata(Some(metadata))
        .set_max_row_group_row_count(Some(group_size))
        .build();
    let mut writer = ArrowWriter::try_new(File::create(&path)?, schema.clone(), Some(properties))?;
    let rows = [
        [2, 0, 1, 1, 0, 0, 1, 1, 1, 1],
        [0, 0, 1, 1, 2, 0, 1, 1, 1, 1],
    ];
    let columns = (0..10)
        .map(|c| Arc::new(Int8Array::from(vec![rows[0][c], rows[1][c]])) as arrow_array::ArrayRef)
        .collect();
    writer.write(&RecordBatch::try_new(schema, columns)?)?;
    writer.close()?;
    seal_fixture(&path, &header)?;
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
    assert_eq!(first.block_lengths(), [1, 1]);
    assert_eq!(
        first.block_symmetries(),
        [BlockSymmetry {
            total_two_j: 0,
            parity: 1
        }; 2]
    );
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
fn reads_only_the_requested_rows_and_checks_their_indexed_symmetry() -> Result<()> {
    let directory = tempfile::tempdir()?;
    // A fully hash-bound pair, but the first block declares odd for an even CSF.
    // The second block is readable even when the two share a row group.
    let blocks =
        "[{two_j=0,parity=-1,start_row=0,length=1},{two_j=0,parity=1,start_row=1,length=1}]";
    for group_size in [1, 2] {
        let (path, header) = fixture_with_index(directory.path(), blocks, group_size, |_| {})?;
        let mut reader = DescriptorReader::open(&path, &header, 16)?;
        let mut values = Vec::new();
        for _ in 0..2 {
            reader.select_block(1)?;
            assert_eq!(reader.read_batch(&mut values)?, 1);
            assert_eq!(values, [0, 0, 1, 1, 2, 0, 1, 1, 1, 1]);
            assert_eq!(reader.read_batch(&mut values)?, 0);
        }
        reader.select_block(0)?;
        let error = reader.read_batch(&mut values).unwrap_err();
        assert!(format!("{error:#}").contains("block 1, CSF 1"));
        assert!(format!("{error:#}").contains("symmetry"));
    }
    Ok(())
}

#[test]
fn rejects_another_valid_parquet_with_the_same_header_body() -> Result<()> {
    let directory = tempfile::tempdir()?;
    let (_, original_header) = fixture(directory.path(), |_| {})?;
    let original_bytes = std::fs::read(original_header)?;
    let (path, header) = fixture(directory.path(), |kv| {
        kv.push(KeyValue::new("extra".into(), Some("changed".into())))
    })?;
    std::fs::write(&header, &original_bytes)?;
    assert!(
        DescriptorReader::open(&path, &header, 1)
            .err()
            .unwrap()
            .to_string()
            .contains("descriptor_sha256")
    );
    let (_, legacy_body) =
        original_bytes.split_at(original_bytes.iter().position(|&b| b == b'\n').unwrap() + 1);
    std::fs::write(&header, legacy_body)?;
    assert!(DescriptorReader::open(&path, &header, 1).is_err());
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
    use crate::descriptor_ffi::{
        rcsfs_descriptor_abi_version, rcsfs_descriptor_close, rcsfs_descriptor_format_version,
        rcsfs_descriptor_info, rcsfs_descriptor_open, rcsfs_descriptor_read,
        rcsfs_descriptor_select_block, rcsfs_descriptor_tables,
    };
    use std::ffi::{CStr, CString, c_char, c_void};
    let directory = tempfile::tempdir()?;
    let (path, header) = fixture(directory.path(), |_| {})?;
    let path = CString::new(path.to_str().unwrap())?;
    let header = CString::new(header.to_str().unwrap())?;
    let mut handles = [std::ptr::null_mut::<c_void>(); 2];
    let mut error = [0 as c_char; 512];
    // SAFETY: The CString paths and stack buffers stay live through each call
    // and each advertised capacity fits the corresponding buffer. Handle
    // slots start null and successful opens are used serially before close.
    // The repeated close uses an initialized null slot, and info checks a null handle.
    unsafe {
        assert_eq!(rcsfs_descriptor_abi_version(), 3);
        assert_eq!(rcsfs_descriptor_format_version(), 3);
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
            [2, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 2, 0, 1, 1, 1, 1]
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
        // The sequential reader buffered both CSFs. Seeking must discard its unread tail.
        for block in [0, 1, 0] {
            assert_eq!(
                rcsfs_descriptor_select_block(handles[1], block, error.as_mut_ptr(), 512),
                0
            );
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
            assert_eq!(other, rows[10 * block..10 * (block + 1)]);
        }
        other.fill(-99);
        for _ in 0..2 {
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
            let error_bytes = error.map(|byte| byte as u8);
            assert_eq!(
                CStr::from_bytes_until_nul(&error_bytes)?.to_bytes(),
                b"unexpected EOF"
            );
            assert_eq!(other, [-99; 10]);
        }
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
    let text = "[header_info]\nheader_lines = ['Core subshells:', '1s', 'Peel subshells:', '', 'CSF(s):']\nelectron_count = 2\n[block_info]\nblock_lengths = [1]\nblocks=[{two_j=0,parity=1,start_row=0,length=1}]\n";
    std::fs::write(&header, text)?;
    let layout = DescriptorLayout { subshell_count: 0 };
    let schema = Arc::new(layout.schema());
    let kv = fixture_metadata(layout, &[], &header)?;
    let properties = WriterProperties::builder()
        .set_key_value_metadata(Some(kv))
        .build();
    let mut writer = ArrowWriter::try_new(File::create(&path)?, schema.clone(), Some(properties))?;
    writer.write(&RecordBatch::try_new(
        schema,
        vec![
            Arc::new(Int8Array::from(vec![1])),
            Arc::new(Int8Array::from(vec![1])),
        ],
    )?)?;
    writer.close()?;
    seal_fixture(&path, &header)?;
    let mut reader = DescriptorReader::open(&path, &header, 1)?;
    let mut values = Vec::new();
    assert_eq!(reader.read_batch(&mut values)?, 1);
    assert_eq!(values, [1, 1]);
    assert_eq!(reader.read_batch(&mut values)?, 0);
    assert!(values.is_empty());
    std::fs::write(&header, format!("{text}\n"))?;
    assert!(DescriptorReader::open(&path, &header, 1).is_err());
    assert!(DescriptorReader::open(&path, &header, 0).is_err());
    Ok(())
}
