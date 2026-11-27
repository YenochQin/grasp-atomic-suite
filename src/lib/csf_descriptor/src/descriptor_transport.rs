//! On-disk layout of the Suite's V3 CSF descriptors.
use anyhow::{Result, ensure};
use arrow_schema::{DataType, Field, Schema};

/// On-disk version tag accepted by this reader.
pub const DESCRIPTOR_VERSION: u8 = 3;
/// V3 channels per peel subshell, ordered by occupation, seniority, 2J+1 and 2K+1.
pub const CHANNELS_PER_SUBSHELL: usize = 4;
/// Largest doubled angular momentum encodable as `2J+1` in a signed byte.
pub const MAX_TWO_J: i32 = 126;

/// The V3 row shape for a fixed number of peel subshells.
#[derive(Clone, Copy, Debug)]
pub struct DescriptorLayout {
    /// Number of peel subshells in column order. The closed core has no channels.
    pub subshell_count: usize,
}

impl DescriptorLayout {
    /// Total channels per row, including the two global channels.
    #[must_use]
    pub fn descriptor_row_channel_count(self) -> usize {
        CHANNELS_PER_SUBSHELL * self.subshell_count + 2
    }

    /// Column names in the order of the subshell feature channels.
    #[must_use]
    pub fn feature_column_names(self) -> Vec<String> {
        (0..self.subshell_count)
            .flat_map(|i| {
                [
                    "occupation",
                    "seniority",
                    "two_j_plus_one",
                    "cumulative_two_k_plus_one",
                ]
                .map(|channel| format!("sub{i}_{channel}"))
            })
            .collect()
    }

    /// Names of the final two columns, after all subshell channels.
    #[must_use]
    pub fn global_column_names(self) -> Vec<String> {
        vec!["total_two_j_plus_one".into(), "parity".into()]
    }

    /// Arrow schema for the ordered non-null Int8 descriptor columns.
    #[must_use]
    pub fn schema(self) -> Schema {
        Schema::new(
            self.feature_column_names()
                .into_iter()
                .chain(self.global_column_names())
                .map(|name| Field::new(name, DataType::Int8, false))
                .collect::<Vec<_>>(),
        )
    }

    /// Verify field names, order, types and nullability against the V3 layout.
    pub fn validate_schema(self, actual: &Schema) -> Result<()> {
        ensure!(
            actual.fields() == self.schema().fields(),
            "descriptor schema must contain the ordered non-null Int8 V3 columns"
        );
        Ok(())
    }
}
