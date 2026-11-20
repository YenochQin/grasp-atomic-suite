//! Read V2 descriptor batches without restoring text or invoking Python.
//! Buffers contain original channels; physical validation belongs to the consumer.
use anyhow::{Context, Result, ensure};
use arrow_array::{Array, Int8Array};
use parquet::arrow::arrow_reader::{ParquetRecordBatchReader, ParquetRecordBatchReaderBuilder};
use sha2::{Digest, Sha256};
use std::{fs::File, path::Path};

use crate::descriptor_schema::{DescriptorLayout, DescriptorVersion};

pub struct DescriptorReader {
    reader: ParquetRecordBatchReader,
    pub header_lines: [String; 5],
    pub block_lengths: Vec<u64>,
    pub zero_order_counts: Option<Vec<i64>>,
    pub layout: DescriptorLayout,
    remaining: u64,
    failed: bool,
    pub core_count: usize,
    pub orbitals: Vec<(i32, i32)>,
    pub block_symmetries: Vec<(i32, i32)>,
    pub electron_count: i32,
    block_index: usize,
    block_rows: u64,
}

fn orbital(label: &str) -> Result<(i32, i32)> {
    ensure!(label.is_ascii(), "non-ASCII orbital label");
    let suffix_start = label
        .bytes()
        .position(|b| b.is_ascii_alphabetic())
        .unwrap_or(label.len());
    let (n, suffix) = label.split_at(suffix_start);
    ensure!(
        !n.is_empty() && !suffix.is_empty(),
        "invalid subshell label {label}"
    );
    ensure!(suffix.len() <= 2, "invalid subshell suffix in {label}");
    let l = "spdfghiklmno"
        .find(&suffix[..1])
        .with_context(|| format!("unsupported orbital letter in {label}"))? as i32;
    let n: i32 = n
        .parse()
        .with_context(|| format!("invalid principal quantum number in {label}"))?;
    ensure!(n > l, "principal quantum number must exceed l in {label}");
    let kappa = match suffix.as_bytes().get(1) {
        None => -l - 1,
        Some(b'-') if l > 0 => l,
        _ => anyhow::bail!("invalid relativistic suffix in {label}"),
    };
    Ok((n, kappa))
}

fn read_header(header: &Path) -> Result<(Vec<u8>, toml::Value)> {
    let read = || -> Result<_> {
        let bytes = std::fs::read(header)?;
        let document = toml::from_str(std::str::from_utf8(&bytes)?)?;
        Ok((bytes, document))
    };
    read().with_context(|| format!("cannot parse TOML header {}", header.display()))
}

fn header_info(document: &toml::Value) -> Result<&toml::Table> {
    document
        .get("header_info")
        .and_then(toml::Value::as_table)
        .context("missing or invalid [header_info] table")
}

fn declared_electrons(info: &toml::Table) -> Result<Option<i32>> {
    info.get("electron_count")
        .map(|value| {
            let number = value.as_integer().context(
                "header_info.electron_count must be an integer including core electrons",
            )?;
            i32::try_from(number)
                .ok()
                .filter(|&n| n > 0)
                .context("electron_count must be positive and fit a suite integer")
        })
        .transpose()
}

fn header_lines(info: &toml::Table) -> Result<[&str; 5]> {
    let lines = info
        .get("header_lines")
        .and_then(toml::Value::as_array)
        .context("header_info.header_lines must be an array")?;
    ensure!(
        lines.len() == 5,
        "header_info.header_lines must contain five strings"
    );
    let mut result = [""; 5];
    for (slot, value) in result.iter_mut().zip(lines) {
        let line = value
            .as_str()
            .context("header_info.header_lines must contain strings only")?;
        ensure!(
            !line.contains(['\n', '\r', '\0']),
            "header_lines entries must not contain newline, carriage return or NUL"
        );
        *slot = line;
    }
    Ok(result)
}

const HEADER_LABELS: [(usize, &str); 3] = [
    (0, "Core subshells:"),
    (2, "Peel subshells:"),
    (4, "CSF(s):"),
];

fn header_orbitals(core: &[&str], peel: &[&str]) -> Result<Vec<(i32, i32)>> {
    let orbitals = core
        .iter()
        .chain(peel)
        .map(|s| orbital(s))
        .collect::<Result<Vec<_>>>()?;
    ensure!(!orbitals.is_empty(), "empty orbital table");
    let unique: std::collections::HashSet<_> = orbitals.iter().copied().collect();
    ensure!(
        unique.len() == orbitals.len(),
        "duplicate orbital across core/peel lists"
    );
    Ok(orbitals)
}

#[cfg(feature = "descriptor-ffi")]
pub(crate) struct RadialHeader {
    pub orbitals: Vec<(i32, i32)>,
    pub core_count: usize,
    pub electron_count: i32,
}

#[cfg(feature = "descriptor-ffi")]
impl RadialHeader {
    pub fn read(path: &Path) -> Result<Self> {
        let (_, document) = read_header(path)?;
        let info = header_info(&document)?;
        let electron_count = declared_electrons(info)?
            .context("header_info.electron_count is required; regenerate old headers with rCSFs")?;
        let lines = header_lines(info)?;
        for (i, expected) in HEADER_LABELS {
            ensure!(
                lines[i].trim_end_matches(' ') == expected,
                "invalid header label"
            );
        }
        // Preserve the separators accepted by the radial input, including commas.
        let [core, peel] = [1, 3].map(|i| {
            lines[i]
                .split([' ', ',', '\t'])
                .filter(|s| !s.is_empty())
                .collect::<Vec<_>>()
        });
        for label in core.iter().chain(&peel) {
            ensure!(
                label.as_bytes()[0].is_ascii_digit(),
                "invalid subshell label {label}"
            );
        }
        Ok(Self {
            orbitals: header_orbitals(&core, &peel)?,
            core_count: core.len(),
            electron_count,
        })
    }
}

impl DescriptorReader {
    /// Open with a bounded output batch size. Header hashing and parsing use
    /// the same byte snapshot. No application global state is modified.
    pub fn open(path: &Path, header: &Path, batch_size: usize) -> Result<Self> {
        ensure!(
            batch_size > 0 && batch_size <= 8192,
            "batch size must be in 1..=8192"
        );
        let (bytes, document) = read_header(header)?;
        let digest = format!("{:x}", Sha256::digest(&bytes));
        let info = header_info(&document)?;
        let header_lines = header_lines(info)?;
        for (i, expected) in HEADER_LABELS {
            ensure!(header_lines[i].trim() == expected, "invalid header label");
        }
        let peel: Vec<_> = header_lines[3].split_whitespace().collect();
        let core: Vec<_> = header_lines[1].split_whitespace().collect();
        let orbitals = header_orbitals(&core, &peel)?;
        let core_count = core.len();
        let layout = DescriptorLayout::new(DescriptorVersion::V2, peel.len());
        let lengths = document
            .get("block_info")
            .and_then(|v| v.get("block_lengths"))
            .and_then(toml::Value::as_array)
            .context("missing block lengths")?;
        let block_lengths = lengths
            .iter()
            .map(|v| {
                let n = v.as_integer().context("block length is not an integer")?;
                ensure!(n > 0, "empty or negative block");
                Ok(n as u64)
            })
            .collect::<Result<Vec<_>>>()?;
        ensure!(!block_lengths.is_empty(), "no blocks");
        let zero_order_counts = document
            .get("block_info")
            .and_then(|v| v.get("zero_order_counts"))
            .map(|v| {
                let counts = v.as_array().context("block_info.zero_order_counts must be an array")?;
                ensure!(
                    counts.len() == block_lengths.len(),
                    "block_info.zero_order_counts length does not match block count"
                );
                counts.iter().zip(&block_lengths).enumerate().map(|(i, (value, &length))| {
                    value.as_integer()
                        .filter(|&n| n > 0 && (n as u64) <= length)
                        .with_context(|| format!("invalid zero-order CSF count for block {} (expected 1..={length})", i + 1))
                }).collect::<Result<Vec<_>>>()
            })
            .transpose()?;
        let total = block_lengths.iter().try_fold(0u64, |sum, n| {
            sum.checked_add(*n).context("row count overflow")
        })?;
        for (section, key, expected) in [
            ("block_info", "block_count", block_lengths.len() as u64),
            ("conversion_stats", "csf_count", total),
        ] {
            if let Some(value) = document.get(section).and_then(|v| v.get(key)) {
                ensure!(
                    value.as_integer().and_then(|n| u64::try_from(n).ok()) == Some(expected),
                    "inconsistent {key}"
                );
            }
        }
        let file = File::open(path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file.try_clone()?)?;
        let metadata = builder.metadata().file_metadata();
        ensure!(
            u64::try_from(metadata.num_rows()).ok() == Some(total),
            "row count mismatch"
        );
        let kv = metadata
            .key_value_metadata()
            .context("missing descriptor metadata")?;
        let get = |key: &str| -> Result<&str> {
            let mut entries = kv.iter().filter(|v| v.key == key);
            let value = entries
                .next()
                .and_then(|v| v.value.as_deref())
                .context(format!("missing {key}"))?;
            ensure!(entries.next().is_none(), "duplicate {key}");
            Ok(value)
        };
        for (key, expected) in [
            ("descriptor_version", "2"),
            ("channels_per_subshell", "4"),
            ("missing_sentinel", "-1"),
            ("normalized", "false"),
            ("source_header_sha256", digest.as_str()),
        ] {
            ensure!(get(key)? == expected, "invalid {key}");
        }
        ensure!(
            get("subshell_count")?.parse::<usize>()? == peel.len(),
            "peel count mismatch"
        );
        crate::descriptor_schema::validate_seniority_encoding(kv)?;
        crate::descriptor_schema::validate_descriptor_storage(kv)?;
        crate::descriptor_schema::validate_output_schema(builder.schema(), layout)?;
        for (key, expected) in [
            (
                "peel_subshells",
                peel.iter().map(|s| s.to_string()).collect::<Vec<_>>(),
            ),
            ("feature_columns", layout.feature_column_names()),
            ("global_columns", layout.global_column_names()),
        ] {
            ensure!(
                serde_json::from_str::<Vec<String>>(get(key)?)? == expected,
                "invalid {key}"
            );
        }
        let declared_electrons = declared_electrons(info)?;
        let reader = builder.with_batch_size(batch_size).build()?;
        let mut source = Self {
            reader,
            header_lines: header_lines.map(str::to_owned),
            block_lengths,
            zero_order_counts,
            layout,
            remaining: total,
            failed: false,
            core_count,
            orbitals,
            block_symmetries: Vec::new(),
            electron_count: declared_electrons.unwrap_or(0),
            block_index: 0,
            block_rows: 0,
        };
        // Scan once to establish each block's symmetry and the electron count before Fortran
        // installs the basis header. Reopen the stream for delivery with the tables in place.
        let mut values = Vec::new();
        while source.read_batch(&mut values)? != 0 {}
        ensure!(
            source.remaining == 0
                && source.block_index == source.block_lengths.len()
                && source.block_rows == 0
                && source.block_symmetries.len() == source.block_lengths.len(),
            "incomplete descriptor block table"
        );
        source.reader = ParquetRecordBatchReaderBuilder::try_new(file)?
            .with_batch_size(batch_size)
            .build()?;
        source.remaining = total;
        source.block_index = 0;
        source.block_rows = 0;
        Ok(source)
    }

    /// Returns rows copied into caller-owned channel-major Fortran storage
    /// (`channel + row_width * row`). An error invalidates this reader.
    pub fn read_batch(&mut self, output: &mut Vec<i8>) -> Result<usize> {
        output.clear();
        ensure!(!self.failed, "reader failed; close and reopen it");
        self.failed = true;
        let batch = loop {
            let Some(batch) = self.reader.next() else {
                ensure!(self.remaining == 0, "premature end of descriptors");
                self.failed = false;
                return Ok(0);
            };
            let batch = batch?;
            if batch.num_rows() > 0 {
                break batch;
            }
        };
        ensure!(
            batch.num_columns() == self.layout.row_len(),
            "column count mismatch"
        );
        let columns = batch
            .columns()
            .iter()
            .map(|column| {
                ensure!(column.null_count() == 0, "null descriptor channel");
                column
                    .as_any()
                    .downcast_ref::<Int8Array>()
                    .context("non-Int8 descriptor channel")
            })
            .collect::<Result<Vec<_>>>()?;
        let count = batch.num_rows();
        self.remaining = self
            .remaining
            .checked_sub(count as u64)
            .context("too many rows")?;
        let size = count
            .checked_mul(self.layout.row_len())
            .context("batch size overflow")?;
        output.try_reserve(size)?;
        // Transpose the Arrow column layout into contiguous descriptor rows for Fortran.
        for row in 0..count {
            for column in &columns {
                output.push(column.value(row));
            }
        }
        for row in output.chunks_exact(self.layout.row_len()) {
            ensure!(
                self.block_index < self.block_lengths.len(),
                "row outside blocks"
            );
            let symmetry = (i32::from(row[row.len() - 2]), i32::from(row[row.len() - 1]));
            ensure!(
                (0..=126).contains(&symmetry.0) && symmetry.1.abs() == 1,
                "invalid symmetry"
            );
            if self.block_symmetries.len() == self.block_index {
                self.block_symmetries.push(symmetry);
            }
            ensure!(
                self.block_symmetries[self.block_index] == symmetry,
                "mixed block symmetry"
            );
            let mut electrons: i32 = self.orbitals[..self.core_count]
                .iter()
                .map(|(_, k)| 2 * k.abs())
                .sum();
            // The core is filled, so it contributes electrons and an even parity exponent.
            let mut parity_sum = 0;
            for (shell, &(_, k)) in self.orbitals[self.core_count..].iter().enumerate() {
                let q = i32::from(row[4 * shell]);
                ensure!(
                    (0..=2 * k.abs()).contains(&q),
                    "occupation outside capacity"
                );
                ensure!(row[4 * shell + 2] >= 0, "negative seniority");
                ensure!(
                    (-1..=126).contains(&row[4 * shell + 1])
                        && (-1..=126).contains(&row[4 * shell + 3]),
                    "invalid J/K"
                );
                if q == 0 {
                    ensure!(
                        row[4 * shell..4 * shell + 4] == [0, -1, 0, -1],
                        "noncanonical empty shell"
                    );
                }
                electrons += q;
                parity_sum += q * if k > 0 { k } else { -k - 1 };
            }
            ensure!(electrons > 0, "empty electron basis");
            if self.electron_count == 0 {
                self.electron_count = electrons;
            }
            ensure!(electrons == self.electron_count, "electron count mismatch");
            ensure!(
                symmetry.1 == if parity_sum % 2 == 0 { 1 } else { -1 },
                "incorrect parity"
            );
            if self.layout.subshell_count() == 0 {
                ensure!(
                    self.block_lengths == [1] && symmetry == (0, 1),
                    "invalid closed-core basis"
                );
            }
            self.block_rows += 1;
            if self.block_rows == self.block_lengths[self.block_index] {
                self.block_index += 1;
                self.block_rows = 0;
            }
        }
        self.failed = false;
        Ok(count)
    }
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
