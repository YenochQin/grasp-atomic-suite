//! Validate and stream native V3 descriptor rows.
use anyhow::{Context, Result, ensure};
use arrow_array::{Array, Int8Array};
use parquet::arrow::arrow_reader::{
    ArrowReaderMetadata, ArrowReaderOptions, ParquetRecordBatchReader,
    ParquetRecordBatchReaderBuilder,
};
use sha2::{Digest, Sha256};
use std::{fs::File, io::Seek, path::Path, sync::LazyLock};

use crate::descriptor_schema::{
    CHANNELS_PER_SUBSHELL, DESCRIPTOR_VERSION, DescriptorLayout, MAX_TWO_J,
};

const ORBITAL_LETTERS: &str = "spdfghiklmn";

pub struct DescriptorReader {
    reader: Option<ParquetRecordBatchReader>,
    file: File,
    metadata: ArrowReaderMetadata,
    batch_size: usize,
    // Cumulative logical row offsets, including the end sentinel.
    row_group_offsets: Vec<u64>,
    block_starts: Vec<u64>,
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
    // Physical states may be readable by consumers that do not use the spin-angular library.
    require_spin_angular: bool,
    block_index: usize,
    block_rows: u64,
}

fn parse_orbital_label(label: &str) -> Result<(i32, i32)> {
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
    let l = ORBITAL_LETTERS
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

fn declared_electrons(info: &toml::Table) -> Result<i32> {
    let number = info
        .get("electron_count")
        .context("header_info.electron_count is required; regenerate old headers with rCSFs")?
        .as_integer()
        .context("header_info.electron_count must be an integer including core electrons")?;
    i32::try_from(number)
        .ok()
        .filter(|&n| n > 0)
        .context("electron_count must be positive and fit a suite integer")
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
        .map(|s| parse_orbital_label(s))
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
        let electron_count = declared_electrons(info)?;
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
        let expected = document
            .get("descriptor_sha256")
            .and_then(toml::Value::as_str)
            .context("missing descriptor_sha256; regenerate the indexed descriptor pair")?;
        ensure!(
            expected.len() == 64
                && expected
                    .bytes()
                    .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b)),
            "invalid descriptor_sha256"
        );
        let prefix = format!("descriptor_sha256 = \"{expected}\"\n");
        let body = bytes
            .strip_prefix(prefix.as_bytes())
            .context("invalid descriptor_sha256 header line")?;
        let digest = format!("{:x}", Sha256::digest(body));
        let info = header_info(&document)?;
        let header_lines = header_lines(info)?;
        for (i, expected) in HEADER_LABELS {
            ensure!(header_lines[i].trim() == expected, "invalid header label");
        }
        let peel: Vec<_> = header_lines[3].split_whitespace().collect();
        let core: Vec<_> = header_lines[1].split_whitespace().collect();
        let orbitals = header_orbitals(&core, &peel)?;
        let core_count = core.len();
        let layout = DescriptorLayout {
            subshell_count: peel.len(),
        };
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
                        .filter(|&n| n >= 0 && (n as u64) <= length)
                        .with_context(|| format!("invalid zero-order CSF count for block {} (expected 0..={length})", i + 1))
                }).collect::<Result<Vec<_>>>()
            })
            .transpose()?;
        let blocks = document["block_info"]
            .get("blocks")
            .and_then(toml::Value::as_array)
            .context("missing block_info.blocks index")?;
        ensure!(
            blocks.len() == block_lengths.len(),
            "block index count mismatch"
        );
        let mut block_starts = Vec::with_capacity(blocks.len());
        let mut block_symmetries = Vec::with_capacity(blocks.len());
        let mut total = 0u64;
        for (block, &length) in blocks.iter().zip(&block_lengths) {
            let integer = |name: &str| {
                block
                    .get(name)
                    .and_then(toml::Value::as_integer)
                    .with_context(|| format!("missing or invalid block {name}"))
            };
            let start = u64::try_from(integer("start_row")?)?;
            ensure!(
                start == total && u64::try_from(integer("length")?)? == length,
                "invalid block index range"
            );
            let two_j = i32::try_from(integer("two_j")?)?;
            let parity = i32::try_from(integer("parity")?)?;
            ensure!(
                (0..=MAX_TWO_J).contains(&two_j) && matches!(parity, -1 | 1),
                "invalid indexed symmetry"
            );
            block_starts.push(start);
            block_symmetries.push((two_j, parity));
            total = total.checked_add(length).context("row count overflow")?;
        }
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
        let mut file = File::open(path)?;
        // Derive the logical types from the Parquet schema, not an embedded Arrow schema.
        let arrow_metadata = ArrowReaderMetadata::load(
            &file,
            ArrowReaderOptions::new().with_skip_arrow_metadata(true),
        )?;
        let builder = ParquetRecordBatchReaderBuilder::new_with_metadata(
            file.try_clone()?,
            arrow_metadata.clone(),
        );
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
        let version = get("descriptor_version")?;
        ensure!(
            version == DESCRIPTOR_VERSION.to_string(),
            "unsupported descriptor_version={version}; regenerate V3 descriptors from the original CSF data"
        );
        ensure!(
            !kv.iter().any(|entry| entry.key == "missing_sentinel"),
            "V3 metadata must omit missing_sentinel"
        );
        for (key, expected) in [
            ("channels_per_subshell", "4"),
            ("descriptor_dtype", "int8"),
            ("max_two_j", "126"),
            ("seniority_encoding", "physical"),
            ("normalized", "false"),
            ("source_header_sha256", digest.as_str()),
        ] {
            ensure!(get(key)? == expected, "invalid {key}");
        }
        ensure!(
            get("subshell_count")?.parse::<usize>()? == peel.len(),
            "peel count mismatch"
        );
        layout.validate_schema(builder.schema())?;
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
        let electron_count = declared_electrons(info)?;
        file.rewind()?;
        let mut file_digest = Sha256::new();
        std::io::copy(&mut file, &mut file_digest)?;
        ensure!(
            format!("{:x}", file_digest.finalize()) == expected,
            "descriptor_sha256 does not match the final Parquet file"
        );
        file.rewind()?;
        let mut row_group_offsets = vec![0u64];
        for group in builder.metadata().row_groups() {
            let end = row_group_offsets
                .last()
                .unwrap()
                .checked_add(u64::try_from(group.num_rows())?)
                .context("row group count overflow")?;
            row_group_offsets.push(end);
        }
        ensure!(
            row_group_offsets.last() == Some(&total),
            "row group count mismatch"
        );
        Ok(Self {
            reader: None,
            file,
            metadata: arrow_metadata,
            batch_size,
            row_group_offsets,
            block_starts,
            header_lines: header_lines.map(str::to_owned),
            block_lengths,
            zero_order_counts,
            layout,
            remaining: total,
            failed: false,
            core_count,
            orbitals,
            block_symmetries,
            electron_count,
            block_index: 0,
            require_spin_angular: false,
            block_rows: 0,
        })
    }

    /// Whole-basis angular programs need this decision before creating MCP files.
    /// Check only high-j shells whose occupation statistics overlap the unsupported range.
    pub fn check_spin_angular(&mut self) -> Result<()> {
        use parquet::arrow::ProjectionMask;
        use parquet::file::statistics::Statistics;

        let metadata = self.metadata.metadata();
        for (shell, &(n, kappa)) in self.orbitals[self.core_count..].iter().enumerate() {
            let capacity = 2 * kappa.abs();
            if capacity < 10 {
                continue;
            }
            let occupation_column = shell * CHANNELS_PER_SUBSHELL;
            for (group_id, group) in metadata.row_groups().iter().enumerate() {
                if let Some(Statistics::Int32(stats)) = group.column(occupation_column).statistics()
                    && (stats.max_opt().is_some_and(|&max| max <= 2)
                        || stats.min_opt().is_some_and(|&min| min >= capacity))
                {
                    continue;
                }
                let projection = ProjectionMask::leaves(
                    metadata.file_metadata().schema_descr(),
                    occupation_column..occupation_column + 3,
                );
                let batches = ParquetRecordBatchReaderBuilder::new_with_metadata(
                    self.file.try_clone()?,
                    self.metadata.clone(),
                )
                .with_projection(projection)
                .with_row_groups(vec![group_id])
                .with_batch_size(self.batch_size)
                .build()?;
                let mut start = self.row_group_offsets[group_id];
                for batch in batches {
                    let batch = batch?;
                    let column = |i| {
                        batch
                            .column(i)
                            .as_any()
                            .downcast_ref::<Int8Array>()
                            .context("non-Int8 high-j column")
                    };
                    let occupations = column(0)?;
                    for (i, value) in occupations.iter().enumerate() {
                        let occupation = i32::from(value.context("null high-j occupation")?);
                        if (3..capacity).contains(&occupation) {
                            let ordinal = start + i as u64;
                            let block =
                                self.block_starts.partition_point(|&row| row <= ordinal) - 1;
                            let row = ordinal - self.block_starts[block] + 1;
                            let seniority = i32::from(column(1)?.value(i));
                            let two_j = decode_angular(column(2)?.value(i))?;
                            validate_subshell(kappa.abs(), occupation, seniority, two_j)
                                .with_context(|| {
                                    format!(
                                        "block {}, CSF {row}, orbital n={n}, kappa={kappa}",
                                        block + 1
                                    )
                                })?;
                            anyhow::bail!(
                                "block {}, CSF {row}, orbital n={n}, kappa={kappa}: open j={}/2 shell with occupation {occupation} exceeds the spin-angular term table",
                                block + 1,
                                capacity - 1
                            );
                        }
                    }
                    start += batch.num_rows() as u64;
                }
            }
        }
        self.require_spin_angular = true;
        Ok(())
    }

    pub fn select_block(&mut self, index: usize) -> Result<()> {
        ensure!(!self.failed, "reader failed; close and reopen it");
        let &start = self
            .block_starts
            .get(index)
            .context("block index out of range")?;
        let length = self.block_lengths[index];
        let first = self
            .row_group_offsets
            .partition_point(|&offset| offset <= start)
            - 1;
        let end = self
            .row_group_offsets
            .partition_point(|&offset| offset < start + length);
        let reader = ParquetRecordBatchReaderBuilder::new_with_metadata(
            self.file.try_clone()?,
            self.metadata.clone(),
        )
        .with_row_groups((first..end).collect())
        .with_offset(usize::try_from(start - self.row_group_offsets[first])?)
        .with_limit(usize::try_from(length)?)
        .with_batch_size(self.batch_size)
        .build()?;
        self.reader = Some(reader);
        self.remaining = length;
        self.block_index = index;
        self.block_rows = 0;
        Ok(())
    }

    /// Returns rows copied into caller-owned channel-major Fortran storage
    /// (`channel + row_width * row`). An error invalidates this reader.
    pub fn read_batch(&mut self, output: &mut Vec<i8>) -> Result<usize> {
        output.clear();
        ensure!(!self.failed, "reader failed; close and reopen it");
        self.failed = true;
        // The sequential C/Rust API still starts at row zero if no block was selected.
        if self.reader.is_none() {
            self.reader = Some(
                ParquetRecordBatchReaderBuilder::new_with_metadata(
                    self.file.try_clone()?,
                    self.metadata.clone(),
                )
                .with_batch_size(self.batch_size)
                .build()?,
            );
        }
        let reader = self.reader.as_mut().expect("initialized descriptor reader");
        let batch = loop {
            let Some(batch) = reader.next() else {
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
            batch.num_columns() == self.layout.descriptor_row_channel_count(),
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
            .checked_mul(self.layout.descriptor_row_channel_count())
            .context("batch size overflow")?;
        output.try_reserve(size)?;
        // Transpose the Arrow column layout into contiguous descriptor rows for Fortran.
        for row in 0..count {
            for column in &columns {
                output.push(column.value(row));
            }
        }
        for row in output.chunks_exact(self.layout.descriptor_row_channel_count()) {
            ensure!(
                self.block_index < self.block_lengths.len(),
                "row outside blocks"
            );
            self.validate_row(row).with_context(|| {
                format!(
                    "block {}, CSF {}",
                    self.block_index + 1,
                    self.block_rows + 1
                )
            })?;
            self.block_rows += 1;
            if self.block_rows == self.block_lengths[self.block_index] {
                self.block_index += 1;
                self.block_rows = 0;
            }
        }
        self.failed = false;
        Ok(count)
    }

    fn validate_row(&mut self, row: &[i8]) -> Result<()> {
        let (peel, globals) = row.split_at(row.len() - 2);
        let total_two_j = decode_angular(globals[0])?;
        let parity = i32::from(globals[1]);
        let symmetry = (total_two_j, parity);
        ensure!(matches!(parity, 1 | -1), "invalid parity code");
        ensure!(
            self.block_symmetries[self.block_index] == symmetry,
            "mixed block symmetry"
        );
        if peel.is_empty() {
            ensure!(
                self.block_lengths == [1] && symmetry == (0, 1),
                "invalid closed-core basis"
            );
        }
        let mut electrons: i32 = self.orbitals[..self.core_count]
            .iter()
            .map(|(_, k)| 2 * k.abs())
            .sum();
        let mut parity_exponent = 0;
        let mut previous_two_k = 0;
        let (subshell_rows, remainder) = peel.as_chunks::<CHANNELS_PER_SUBSHELL>();
        debug_assert!(remainder.is_empty());
        for (channels, &(n, kappa)) in subshell_rows.iter().zip(&self.orbitals[self.core_count..]) {
            let (occupation, seniority) = (i32::from(channels[0]), i32::from(channels[1]));
            let validate_shell = || -> Result<_> {
                let two_j = decode_angular(channels[2])?;
                let two_k = decode_angular(channels[3])?;
                validate_subshell(kappa.abs(), occupation, seniority, two_j)?;
                let maximum = previous_two_k + two_j;
                ensure!(
                    ((previous_two_k - two_j).abs()..=maximum).contains(&two_k)
                        && (maximum - two_k) % 2 == 0,
                    "cumulative coupling violates angular triangle or step"
                );
                Ok(two_k)
            };
            previous_two_k =
                validate_shell().with_context(|| format!("orbital n={n}, kappa={kappa}"))?;
            if kappa.abs() >= 5
                && (3..2 * kappa.abs()).contains(&occupation)
                && self.require_spin_angular
            {
                anyhow::bail!(
                    "block {}, CSF {}, orbital n={n}, kappa={kappa}: open j={}/2 shell with occupation {occupation} exceeds the spin-angular term table",
                    self.block_index + 1,
                    self.block_rows + 1,
                    2 * kappa.abs() - 1
                );
            }
            electrons += occupation;
            parity_exponent += occupation * if kappa > 0 { kappa } else { -kappa - 1 };
        }
        ensure!(
            previous_two_k == total_two_j,
            "final cumulative coupling differs from total J"
        );
        ensure!(electrons > 0, "empty electron basis");
        ensure!(electrons == self.electron_count, "electron count mismatch");
        ensure!(
            parity == if parity_exponent % 2 == 0 { 1 } else { -1 },
            "incorrect parity"
        );
        Ok(())
    }
}

fn decode_angular(code: i8) -> Result<i32> {
    let two_j = i32::from(code) - 1;
    ensure!(
        (0..=MAX_TWO_J).contains(&two_j),
        "angular code must be in 1..=127"
    );
    Ok(two_j)
}

fn validate_subshell(
    half_capacity: i32,
    occupation: i32,
    seniority: i32,
    two_j: i32,
) -> Result<()> {
    let capacity = 2 * half_capacity;
    ensure!(
        (0..=capacity).contains(&occupation),
        "occupation outside capacity"
    );
    let population = occupation.min(capacity - occupation);
    ensure!(
        (0..=population).contains(&seniority) && (population - seniority) % 2 == 0,
        "invalid subshell seniority for occupation"
    );
    if population == 0 {
        ensure!(two_j == 0, "empty or full shell requires J=0");
        return Ok(());
    }
    ensure!(
        half_capacity <= 5 || population <= 2,
        "unsupported open-shell occupation"
    );
    let multiplicities = &SUBSHELL_TERMS[(half_capacity - 1) as usize];
    let seniority = seniority as usize;
    let two_j = two_j as usize;
    let angular_multiplicity =
        |particles: usize| multiplicities[particles][two_j] - multiplicities[particles][two_j + 2];
    let paired_multiplicity = seniority.checked_sub(2).map_or(0, angular_multiplicity);
    match angular_multiplicity(seniority) - paired_multiplicity {
        1 => Ok(()),
        0 => anyhow::bail!("invalid subshell (occupation, seniority, J)"),
        _ => anyhow::bail!("unsupported subshell state requires an additional w label"),
    }
}

// Count determinants at each magnetic projection for the particle/hole populations
// supported by Suite's subshell term table. The difference at M=J and M=J+1
// counts J terms. Subtracting the contribution of v-2 removes the paired terms.
const PROJECTION_COUNT: usize = MAX_TWO_J as usize + 3;
const PROJECTION_ORIGIN: usize = PROJECTION_COUNT - 1;
type ProjectionCounts = [i32; PROJECTION_COUNT];
static SUBSHELL_TERMS: LazyLock<Vec<Vec<ProjectionCounts>>> = LazyLock::new(|| {
    (1..=ORBITAL_LETTERS.len() as isize)
        .map(|half_capacity| {
            let single_two_j = 2 * half_capacity - 1;
            let max_population = if half_capacity <= 5 { half_capacity } else { 2 };
            let mut counts = vec![[0; 2 * PROJECTION_ORIGIN + 1]; max_population as usize + 1];
            counts[0][PROJECTION_ORIGIN] = 1;
            for two_m in (-single_two_j..=single_two_j).step_by(2) {
                for particles in (1..counts.len()).rev() {
                    for from in 0..counts[particles - 1].len() {
                        let to = from as isize + two_m;
                        if (0..counts[particles].len() as isize).contains(&to) {
                            counts[particles][to as usize] += counts[particles - 1][from];
                        }
                    }
                }
            }
            counts
                .iter()
                .map(|projections| {
                    projections[PROJECTION_ORIGIN..PROJECTION_ORIGIN + PROJECTION_COUNT]
                        .try_into()
                        .unwrap()
                })
                .collect()
        })
        .collect()
});

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
