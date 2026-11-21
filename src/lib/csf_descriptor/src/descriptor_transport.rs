//! On-disk layout of the Suite's V3 CSF descriptors.
use anyhow::{Result, ensure};
use arrow_schema::{DataType, Field, Schema};

pub const DESCRIPTOR_VERSION: u8 = 3;
pub const CHANNELS_PER_SUBSHELL: usize = 4;
pub const MAX_TWO_J: i32 = 126;

#[derive(Clone, Copy, Debug)]
pub struct DescriptorLayout {
    pub subshell_count: usize,
}

impl DescriptorLayout {
    pub fn descriptor_row_channel_count(self) -> usize {
        CHANNELS_PER_SUBSHELL * self.subshell_count + 2
    }

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

    pub fn global_column_names(self) -> Vec<String> {
        vec!["total_two_j_plus_one".into(), "parity".into()]
    }

    pub fn schema(self) -> Schema {
        Schema::new(
            self.feature_column_names()
                .into_iter()
                .chain(self.global_column_names())
                .map(|name| Field::new(name, DataType::Int8, false))
                .collect::<Vec<_>>(),
        )
    }

    pub fn validate_schema(self, actual: &Schema) -> Result<()> {
        ensure!(
            actual.fields() == self.schema().fields(),
            "descriptor schema must contain the ordered non-null Int8 V3 columns"
        );
        Ok(())
    }
}
