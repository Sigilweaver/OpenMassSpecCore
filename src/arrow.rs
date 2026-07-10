//! Apache Arrow bridge for [`crate::SpectrumRecord`].
//!
//! This module is gated behind the `arrow` Cargo feature. It exposes:
//!
//! * [`spectrum_record_schema`]: the canonical Arrow [`Schema`] for a
//!   stream of `SpectrumRecord`s.
//! * [`SpectrumBatchBuilder`]: a streaming builder that accumulates rows
//!   and produces a [`RecordBatch`] when finalized.
//!
//! The schema is intentionally flat (one row per spectrum). Peak arrays
//! are stored as `LargeList<Float64>` / `LargeList<Float32>` so that
//! individual spectra can exceed `i32::MAX` peaks without truncation.
//! Precursor fields are inlined as nullable scalar columns; an MS1
//! spectrum has all precursor columns null. This shape is friendly to
//! Polars / DataFusion / DuckDB consumers.
//!
//! The Arrow schema is part of this crate's stable surface for the
//! purposes of consumers that pin `openmassspec-core` directly; any column
//! addition is a minor-version bump, any removal or rename is breaking.

use std::sync::Arc;

use arrow_array::builder::{
    ArrayBuilder, Float32Builder, Float64Builder, Int32Builder, LargeListBuilder, StringBuilder,
    UInt32Builder, UInt8Builder,
};
use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{DataType, Field, Schema, SchemaRef};

use crate::{Activation, Analyzer, MobilityArrayKind, Polarity, ScanMode, SpectrumRecord};

/// Return the canonical [`Schema`] for a `RecordBatch` of spectra.
pub fn spectrum_record_schema() -> SchemaRef {
    let mz_item = Arc::new(Field::new("item", DataType::Float64, false));
    let int_item = Arc::new(Field::new("item", DataType::Float32, false));
    let mob_item = Arc::new(Field::new("item", DataType::Float32, false));
    Arc::new(Schema::new(vec![
        Field::new("index", DataType::UInt32, false),
        Field::new("scan_number", DataType::UInt32, false),
        Field::new("native_id", DataType::Utf8, false),
        Field::new("ms_level", DataType::UInt8, false),
        Field::new("polarity", DataType::Utf8, true),
        Field::new("scan_mode", DataType::Utf8, true),
        Field::new("analyzer", DataType::Utf8, true),
        Field::new("filter", DataType::Utf8, true),
        Field::new("retention_time_sec", DataType::Float64, false),
        Field::new("total_ion_current", DataType::Float64, true),
        Field::new("base_peak_mz", DataType::Float64, true),
        Field::new("base_peak_intensity", DataType::Float64, true),
        Field::new("low_mz", DataType::Float64, true),
        Field::new("high_mz", DataType::Float64, true),
        Field::new("ion_injection_time_ms", DataType::Float64, true),
        Field::new("inv_mobility", DataType::Float64, true),
        Field::new("precursor_target_mz", DataType::Float64, true),
        Field::new("precursor_selected_mz", DataType::Float64, true),
        Field::new("precursor_isolation_width", DataType::Float64, true),
        Field::new("precursor_charge", DataType::Int32, true),
        Field::new("precursor_intensity", DataType::Float64, true),
        Field::new("precursor_collision_energy", DataType::Float64, true),
        Field::new("precursor_ce_is_nce", DataType::UInt8, true),
        Field::new("precursor_native_id", DataType::Utf8, true),
        Field::new("precursor_activation", DataType::Utf8, true),
        Field::new("precursor_analyzer", DataType::Utf8, true),
        Field::new_large_list("mz", mz_item, false),
        Field::new_large_list("intensity", int_item, false),
        Field::new_large_list("inv_mobility_per_peak", mob_item, true),
        Field::new("mobility_array_kind", DataType::Utf8, true),
    ]))
}

fn polarity_str(p: Polarity) -> &'static str {
    match p {
        Polarity::Positive => "positive",
        Polarity::Negative => "negative",
    }
}

fn scan_mode_str(m: ScanMode) -> &'static str {
    match m {
        ScanMode::Profile => "profile",
        ScanMode::Centroid => "centroid",
    }
}

fn analyzer_str(a: Analyzer) -> &'static str {
    match a {
        Analyzer::ITMS => "itms",
        Analyzer::TQMS => "tqms",
        Analyzer::SQMS => "sqms",
        Analyzer::TOFMS => "tof",
        Analyzer::FTMS => "ftms",
        Analyzer::Sector => "sector",
    }
}

fn activation_str(a: Activation) -> &'static str {
    match a {
        Activation::CID => "cid",
        Activation::HCD => "hcd",
        Activation::ETD => "etd",
        Activation::ECD => "ecd",
        Activation::UVPD => "uvpd",
        Activation::PQD => "pqd",
        Activation::PD => "pd",
        Activation::SID => "sid",
        Activation::EThcD => "ethcd",
        Activation::IRMPD => "irmpd",
        Activation::MPID => "mpid",
    }
}

fn mobility_kind_str(k: MobilityArrayKind) -> &'static str {
    match k {
        MobilityArrayKind::InverseReducedVsPerCm2 => "inverse_reduced_k0",
        MobilityArrayKind::DriftTimeMilliseconds => "drift_time_ms",
    }
}

/// Streaming builder that accumulates `SpectrumRecord`s and produces a
/// single Arrow [`RecordBatch`] when finalized.
///
/// All rows in a batch share one `mobility_array_kind` value, recorded
/// once at construction. Push rows with [`push`](Self::push); call
/// [`finish`](Self::finish) to materialize the batch.
pub struct SpectrumBatchBuilder {
    schema: SchemaRef,
    mobility_kind: Option<MobilityArrayKind>,
    index: UInt32Builder,
    scan_number: UInt32Builder,
    native_id: StringBuilder,
    ms_level: UInt8Builder,
    polarity: StringBuilder,
    scan_mode: StringBuilder,
    analyzer: StringBuilder,
    filter: StringBuilder,
    retention_time_sec: Float64Builder,
    total_ion_current: Float64Builder,
    base_peak_mz: Float64Builder,
    base_peak_intensity: Float64Builder,
    low_mz: Float64Builder,
    high_mz: Float64Builder,
    ion_injection_time_ms: Float64Builder,
    inv_mobility: Float64Builder,
    precursor_target_mz: Float64Builder,
    precursor_selected_mz: Float64Builder,
    precursor_isolation_width: Float64Builder,
    precursor_charge: Int32Builder,
    precursor_intensity: Float64Builder,
    precursor_collision_energy: Float64Builder,
    precursor_ce_is_nce: UInt8Builder,
    precursor_native_id: StringBuilder,
    precursor_activation: StringBuilder,
    precursor_analyzer: StringBuilder,
    mz: LargeListBuilder<Float64Builder>,
    intensity: LargeListBuilder<Float32Builder>,
    inv_mobility_per_peak: LargeListBuilder<Float32Builder>,
    mobility_array_kind_col: StringBuilder,
}

impl SpectrumBatchBuilder {
    /// Create a new builder. Pass the `mobility_array_kind` from the
    /// source's [`crate::RunMetadata`] so the resulting Arrow batch
    /// carries the unit/CV interpretation alongside the data.
    pub fn new(mobility_kind: Option<MobilityArrayKind>) -> Self {
        Self {
            schema: spectrum_record_schema(),
            mobility_kind,
            index: UInt32Builder::new(),
            scan_number: UInt32Builder::new(),
            native_id: StringBuilder::new(),
            ms_level: UInt8Builder::new(),
            polarity: StringBuilder::new(),
            scan_mode: StringBuilder::new(),
            analyzer: StringBuilder::new(),
            filter: StringBuilder::new(),
            retention_time_sec: Float64Builder::new(),
            total_ion_current: Float64Builder::new(),
            base_peak_mz: Float64Builder::new(),
            base_peak_intensity: Float64Builder::new(),
            low_mz: Float64Builder::new(),
            high_mz: Float64Builder::new(),
            ion_injection_time_ms: Float64Builder::new(),
            inv_mobility: Float64Builder::new(),
            precursor_target_mz: Float64Builder::new(),
            precursor_selected_mz: Float64Builder::new(),
            precursor_isolation_width: Float64Builder::new(),
            precursor_charge: Int32Builder::new(),
            precursor_intensity: Float64Builder::new(),
            precursor_collision_energy: Float64Builder::new(),
            precursor_ce_is_nce: UInt8Builder::new(),
            precursor_native_id: StringBuilder::new(),
            precursor_activation: StringBuilder::new(),
            precursor_analyzer: StringBuilder::new(),
            mz: LargeListBuilder::new(Float64Builder::new()).with_field(Arc::new(Field::new(
                "item",
                DataType::Float64,
                false,
            ))),
            intensity: LargeListBuilder::new(Float32Builder::new())
                .with_field(Arc::new(Field::new("item", DataType::Float32, false))),
            inv_mobility_per_peak: LargeListBuilder::new(Float32Builder::new())
                .with_field(Arc::new(Field::new("item", DataType::Float32, false))),
            mobility_array_kind_col: StringBuilder::new(),
        }
    }

    /// Schema for the batch produced by this builder.
    pub fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    /// Append one spectrum row.
    pub fn push(&mut self, rec: &SpectrumRecord) {
        self.index.append_value(rec.index as u32);
        self.scan_number.append_value(rec.scan_number);
        self.native_id.append_value(&rec.native_id);
        self.ms_level.append_value(rec.ms_level as u8);
        self.polarity.append_option(rec.polarity.map(polarity_str));
        self.scan_mode
            .append_option(rec.scan_mode.map(scan_mode_str));
        self.analyzer.append_option(rec.analyzer.map(analyzer_str));
        self.filter.append_option(rec.filter.as_deref());
        self.retention_time_sec.append_value(rec.retention_time_sec);
        self.total_ion_current.append_option(rec.total_ion_current);
        self.base_peak_mz.append_option(rec.base_peak_mz);
        self.base_peak_intensity
            .append_option(rec.base_peak_intensity);
        self.low_mz.append_option(rec.low_mz);
        self.high_mz.append_option(rec.high_mz);
        self.ion_injection_time_ms
            .append_option(rec.ion_injection_time_ms);
        self.inv_mobility.append_option(rec.inv_mobility);

        match &rec.precursor {
            Some(p) => {
                self.precursor_target_mz.append_option(p.target_mz);
                self.precursor_selected_mz.append_option(p.selected_mz);
                self.precursor_isolation_width
                    .append_option(p.isolation_width);
                self.precursor_charge.append_option(p.charge);
                self.precursor_intensity.append_option(p.intensity);
                self.precursor_collision_energy
                    .append_option(p.collision_energy);
                self.precursor_ce_is_nce.append_value(u8::from(p.ce_is_nce));
                self.precursor_native_id
                    .append_option(p.precursor_native_id.as_deref());
                self.precursor_activation
                    .append_option(p.activation.map(activation_str));
                self.precursor_analyzer
                    .append_option(p.analyzer.map(analyzer_str));
            }
            None => {
                self.precursor_target_mz.append_null();
                self.precursor_selected_mz.append_null();
                self.precursor_isolation_width.append_null();
                self.precursor_charge.append_null();
                self.precursor_intensity.append_null();
                self.precursor_collision_energy.append_null();
                self.precursor_ce_is_nce.append_null();
                self.precursor_native_id.append_null();
                self.precursor_activation.append_null();
                self.precursor_analyzer.append_null();
            }
        }

        for &v in &rec.mz {
            self.mz.values().append_value(v);
        }
        self.mz.append(true);
        for &v in &rec.intensity {
            self.intensity.values().append_value(v);
        }
        self.intensity.append(true);
        match &rec.inv_mobility_per_peak {
            Some(mob) => {
                for &v in mob {
                    self.inv_mobility_per_peak.values().append_value(v);
                }
                self.inv_mobility_per_peak.append(true);
            }
            None => self.inv_mobility_per_peak.append(false),
        }
        self.mobility_array_kind_col
            .append_option(self.mobility_kind.map(mobility_kind_str));
    }

    /// Number of rows accumulated so far.
    pub fn len(&self) -> usize {
        self.index.len()
    }

    /// `true` if no rows have been pushed.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Materialize the accumulated rows into a [`RecordBatch`].
    pub fn finish(mut self) -> Result<RecordBatch, arrow_schema::ArrowError> {
        let arrays: Vec<ArrayRef> = vec![
            Arc::new(self.index.finish()),
            Arc::new(self.scan_number.finish()),
            Arc::new(self.native_id.finish()),
            Arc::new(self.ms_level.finish()),
            Arc::new(self.polarity.finish()),
            Arc::new(self.scan_mode.finish()),
            Arc::new(self.analyzer.finish()),
            Arc::new(self.filter.finish()),
            Arc::new(self.retention_time_sec.finish()),
            Arc::new(self.total_ion_current.finish()),
            Arc::new(self.base_peak_mz.finish()),
            Arc::new(self.base_peak_intensity.finish()),
            Arc::new(self.low_mz.finish()),
            Arc::new(self.high_mz.finish()),
            Arc::new(self.ion_injection_time_ms.finish()),
            Arc::new(self.inv_mobility.finish()),
            Arc::new(self.precursor_target_mz.finish()),
            Arc::new(self.precursor_selected_mz.finish()),
            Arc::new(self.precursor_isolation_width.finish()),
            Arc::new(self.precursor_charge.finish()),
            Arc::new(self.precursor_intensity.finish()),
            Arc::new(self.precursor_collision_energy.finish()),
            Arc::new(self.precursor_ce_is_nce.finish()),
            Arc::new(self.precursor_native_id.finish()),
            Arc::new(self.precursor_activation.finish()),
            Arc::new(self.precursor_analyzer.finish()),
            Arc::new(self.mz.finish()),
            Arc::new(self.intensity.finish()),
            Arc::new(self.inv_mobility_per_peak.finish()),
            Arc::new(self.mobility_array_kind_col.finish()),
        ];
        RecordBatch::try_new(self.schema.clone(), arrays)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{PrecursorInfo, SpectrumRecord};
    use arrow_array::Array;

    fn rec(index: usize, ms_level: u32, n_peaks: usize, with_mob: bool) -> SpectrumRecord {
        let mz: Vec<f64> = (0..n_peaks).map(|i| 100.0 + i as f64).collect();
        let intensity: Vec<f32> = (0..n_peaks).map(|i| 10.0 + i as f32).collect();
        let mobility = if with_mob {
            Some((0..n_peaks).map(|i| 0.5 + i as f32 * 0.01).collect())
        } else {
            None
        };
        SpectrumRecord {
            index,
            scan_number: (index + 1) as u32,
            native_id: format!("scan={}", index + 1),
            ms_level,
            polarity: Some(Polarity::Positive),
            scan_mode: Some(ScanMode::Centroid),
            analyzer: Some(Analyzer::TOFMS),
            filter: None,
            retention_time_sec: index as f64,
            total_ion_current: Some(intensity.iter().map(|&v| v as f64).sum()),
            base_peak_mz: mz.last().copied(),
            base_peak_intensity: intensity.last().map(|&v| v as f64),
            low_mz: mz.first().copied(),
            high_mz: mz.last().copied(),
            ion_injection_time_ms: None,
            inv_mobility: None,
            precursor: if ms_level >= 2 {
                Some(PrecursorInfo {
                    target_mz: Some(500.0),
                    selected_mz: Some(500.5),
                    isolation_width: Some(2.0),
                    charge: Some(2),
                    ..Default::default()
                })
            } else {
                None
            },
            mz,
            intensity,
            inv_mobility_per_peak: mobility,
        }
    }

    #[test]
    fn schema_round_trip() {
        let mut b = SpectrumBatchBuilder::new(Some(MobilityArrayKind::DriftTimeMilliseconds));
        b.push(&rec(0, 1, 3, true));
        b.push(&rec(1, 2, 4, false));
        let batch = b.finish().unwrap();
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.schema().fields().len(), 30);
        let mz_col = batch
            .column_by_name("mz")
            .unwrap()
            .as_any()
            .downcast_ref::<arrow_array::LargeListArray>()
            .unwrap();
        assert_eq!(mz_col.value_length(0), 3);
        assert_eq!(mz_col.value_length(1), 4);
        let mob_col = batch
            .column_by_name("inv_mobility_per_peak")
            .unwrap()
            .as_any()
            .downcast_ref::<arrow_array::LargeListArray>()
            .unwrap();
        assert!(mob_col.is_valid(0));
        assert!(mob_col.is_null(1));
    }
}
