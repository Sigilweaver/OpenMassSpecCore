//! Apache Arrow bridge for [`openmassspec_core::SpectrumRecord`].
//!
//! It exposes:
//!
//! * [`spectrum_record_schema`]: the canonical Arrow [`Schema`] for a
//!   stream of `SpectrumRecord`s.
//! * [`SpectrumBatchBuilder`]: a streaming builder that accumulates rows
//!   and produces a [`RecordBatch`] when finalized.
//!
//! ```
//! use openmassspec_arrow::SpectrumBatchBuilder;
//! use openmassspec_core::SpectrumRecord;
//!
//! let mut builder = SpectrumBatchBuilder::new(None);
//! builder
//!     .push(&SpectrumRecord::new(0, 1, "scan=1").with_peaks(vec![100.0], vec![5.0]))
//!     .unwrap();
//! let batch = builder.finish().unwrap();
//! assert_eq!(batch.num_rows(), 1);
//! ```
//!
//! Enum columns hold the stable labels from each enum's `as_str` method
//! (e.g. [`Analyzer::as_str`]).
//!
//! The schema is intentionally flat (one row per spectrum). Peak arrays
//! are stored as `LargeList<Float64>` / `LargeList<Float32>` so that
//! individual spectra can exceed `i32::MAX` peaks without truncation.
//! Precursor fields are inlined as nullable scalar columns; an MS1
//! spectrum has all precursor columns null. This shape is friendly to
//! Polars / DataFusion / DuckDB consumers.
//!
//! The schema carries `openmassspec.schema_version` in its metadata
//! ([`SCHEMA_VERSION`]); [`check_schema_version`] rejects a schema written by
//! a newer version. New columns are appended at the end, so a reader that
//! looks columns up by name keeps working.
//!
//! # Arrow versions
//!
//! The Arrow types in this API (`RecordBatch`, `SchemaRef`) come from the
//! arrow-rs release selected by a version-named feature: `arrow-60` (the
//! default) uses arrow 60. A new arrow-rs major is supported by adding a
//! feature (`arrow-61`) in a minor release of this crate, not by changing
//! `openmassspec-core`. This crate has its own version line.
//!
//! # Checks
//!
//! [`SpectrumBatchBuilder::push`] runs [`SpectrumRecord::validate`] (array
//! lengths, non-finite values, unsorted m/z, scan windows) and rejects a
//! native id the builder has already seen, so a malformed record is an
//! error ([`BatchError`]) instead of a row that disagrees with itself.
//! Warnings on [`RunMetadata`](openmassspec_core::RunMetadata) and on
//! chromatograms have no column; only spectrum warnings are carried.

#![cfg(feature = "arrow-60")]

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use arrow_array::builder::{
    ArrayBuilder, Float32Builder, Float64Builder, Int32Builder, LargeListBuilder, MapBuilder,
    MapFieldNames, StringBuilder, StructBuilder, UInt32Builder, UInt8Builder,
};
use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{ArrowError, DataType, Field, Fields, Schema, SchemaRef};

use openmassspec_core::{
    Activation, Analyzer, CollisionEnergyUnit, MobilityArrayKind, Polarity, RecordProblem,
    ScanMode, SpectrumRecord,
};

/// Schema metadata key holding [`SCHEMA_VERSION`].
pub const SCHEMA_VERSION_KEY: &str = "openmassspec.schema_version";

/// Version of the column layout written by this crate. It changes only when
/// a column is removed, renamed or retyped; appended columns do not change
/// it.
pub const SCHEMA_VERSION: u32 = 1;

/// Why a record could not be added to a batch, or a batch could not be
/// built or read.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum BatchError {
    /// The record fails Core's checks (see [`SpectrumRecord::validate`]).
    #[error("invalid record {record_id:?}: {problem}")]
    #[non_exhaustive]
    Record {
        /// Native id of the record.
        record_id: String,
        /// What was wrong.
        problem: RecordProblem,
    },
    /// A value does not fit its column (`index` in `u32`, `ms_level` in `u8`).
    #[error("record {record_id:?}: {column} value {value} does not fit its column")]
    #[non_exhaustive]
    ColumnRange {
        /// Native id of the record.
        record_id: String,
        /// Column name.
        column: &'static str,
        /// The value.
        value: u64,
    },
    /// The schema is missing its version or was written by a newer version.
    #[error("unsupported schema version {found:?} (this crate reads up to {supported})")]
    #[non_exhaustive]
    SchemaVersion {
        /// The version string found, `None` if absent.
        found: Option<String>,
        /// The newest version this crate understands.
        supported: u32,
    },
    /// Arrow itself failed.
    #[error("arrow error: {0}")]
    Arrow(#[from] ArrowError),
}

/// Check the `openmassspec.schema_version` metadata of `schema`.
///
/// # Errors
///
/// [`BatchError::SchemaVersion`] when the key is missing, not a number, or
/// larger than [`SCHEMA_VERSION`].
pub fn check_schema_version(schema: &Schema) -> Result<u32, BatchError> {
    let found = schema.metadata().get(SCHEMA_VERSION_KEY).cloned();
    match found.as_deref().and_then(|v| v.parse::<u32>().ok()) {
        Some(v) if v <= SCHEMA_VERSION => Ok(v),
        _ => Err(BatchError::SchemaVersion {
            found,
            supported: SCHEMA_VERSION,
        }),
    }
}

fn struct_item(fields: Vec<Field>) -> Arc<Field> {
    Arc::new(Field::new(
        "item",
        DataType::Struct(Fields::from(fields)),
        false,
    ))
}

fn energy_fields() -> Vec<Field> {
    vec![
        Field::new("value", DataType::Float64, false),
        Field::new("unit", DataType::Utf8, true),
    ]
}

fn window_fields() -> Vec<Field> {
    vec![
        Field::new("lower", DataType::Float64, false),
        Field::new("upper", DataType::Float64, false),
    ]
}

fn warning_fields() -> Vec<Field> {
    vec![
        Field::new("code", DataType::Utf8, false),
        Field::new("message", DataType::Utf8, false),
    ]
}

fn struct_list_builder(
    fields: Vec<Field>,
    builders: Vec<Box<dyn ArrayBuilder>>,
) -> LargeListBuilder<StructBuilder> {
    LargeListBuilder::new(StructBuilder::new(Fields::from(fields.clone()), builders))
        .with_field(struct_item(fields))
}

// Arrow-spec map field names for the `extra` column.
const EXTRA_ENTRIES_FIELD: &str = "entries";
const EXTRA_KEY_FIELD: &str = "key";
const EXTRA_VALUE_FIELD: &str = "value";

/// Return the canonical [`Schema`] for a `RecordBatch` of spectra.
pub fn spectrum_record_schema() -> SchemaRef {
    let mz_item = Arc::new(Field::new("item", DataType::Float64, false));
    let int_item = Arc::new(Field::new("item", DataType::Float32, false));
    let mob_item = Arc::new(Field::new("item", DataType::Float32, false));
    Arc::new(
        Schema::new(vec![
            Field::new("index", DataType::UInt32, false),
            Field::new("scan_number", DataType::UInt32, false),
            Field::new("native_id", DataType::Utf8, false),
            Field::new("ms_level", DataType::UInt8, false),
            Field::new("polarity", DataType::Utf8, true),
            Field::new("scan_mode", DataType::Utf8, true),
            Field::new("analyzer", DataType::Utf8, true),
            Field::new("filter", DataType::Utf8, true),
            Field::new("retention_time_sec", DataType::Float64, true),
            Field::new("total_ion_current", DataType::Float64, true),
            Field::new("base_peak_mz", DataType::Float64, true),
            Field::new("base_peak_intensity", DataType::Float64, true),
            Field::new("low_mz", DataType::Float64, true),
            Field::new("high_mz", DataType::Float64, true),
            Field::new("ion_injection_time_ms", DataType::Float64, true),
            Field::new("inv_mobility", DataType::Float64, true),
            Field::new("faims_cv", DataType::Float64, true),
            Field::new("precursor_target_mz", DataType::Float64, true),
            Field::new("precursor_selected_mz", DataType::Float64, true),
            Field::new("precursor_isolation_width", DataType::Float64, true),
            Field::new("precursor_charge", DataType::Int32, true),
            Field::new("precursor_intensity", DataType::Float64, true),
            Field::new("precursor_collision_energy", DataType::Float64, true),
            Field::new("precursor_collision_energy_unit", DataType::Utf8, true),
            Field::new("precursor_native_id", DataType::Utf8, true),
            Field::new("precursor_activation", DataType::Utf8, true),
            Field::new("precursor_analyzer", DataType::Utf8, true),
            Field::new("precursor_ccs", DataType::Float64, true),
            Field::new_large_list("mz", mz_item, false),
            Field::new_large_list("intensity", int_item, false),
            Field::new_large_list("inv_mobility_per_peak", mob_item, true),
            Field::new("mobility_array_kind", DataType::Utf8, true),
            Field::new("acquisition_event_id", DataType::UInt32, true),
            Field::new_map(
                "extra",
                EXTRA_ENTRIES_FIELD,
                Field::new(EXTRA_KEY_FIELD, DataType::Utf8, false),
                Field::new(EXTRA_VALUE_FIELD, DataType::Utf8, true),
                false,
                false,
            ),
            // Appended columns (schema version 1): everything below can be
            // ignored by a reader that only knows the columns above.
            Field::new_large_list(
                "precursor_collision_energies",
                struct_item(energy_fields()),
                true,
            ),
            Field::new_large_list("scan_windows", struct_item(window_fields()), false),
            Field::new_large_list("warnings", struct_item(warning_fields()), false),
        ])
        .with_metadata(HashMap::from([(
            SCHEMA_VERSION_KEY.to_string(),
            SCHEMA_VERSION.to_string(),
        )])),
    )
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
    faims_cv: Float64Builder,
    precursor_target_mz: Float64Builder,
    precursor_selected_mz: Float64Builder,
    precursor_isolation_width: Float64Builder,
    precursor_charge: Int32Builder,
    precursor_intensity: Float64Builder,
    precursor_collision_energy: Float64Builder,
    precursor_collision_energy_unit: StringBuilder,
    precursor_native_id: StringBuilder,
    precursor_activation: StringBuilder,
    precursor_analyzer: StringBuilder,
    precursor_ccs: Float64Builder,
    mz: LargeListBuilder<Float64Builder>,
    intensity: LargeListBuilder<Float32Builder>,
    inv_mobility_per_peak: LargeListBuilder<Float32Builder>,
    mobility_array_kind_col: StringBuilder,
    acquisition_event_id: UInt32Builder,
    extra: MapBuilder<StringBuilder, StringBuilder>,
    precursor_collision_energies: LargeListBuilder<StructBuilder>,
    scan_windows: LargeListBuilder<StructBuilder>,
    warnings: LargeListBuilder<StructBuilder>,
    /// Native ids already pushed, kept across [`finish_batch`](Self::finish_batch).
    seen_ids: HashSet<String>,
}

impl SpectrumBatchBuilder {
    /// Create a new builder. Pass the `mobility_array_kind` from the
    /// source's [`RunMetadata`](openmassspec_core::RunMetadata) so the resulting Arrow batch
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
            faims_cv: Float64Builder::new(),
            precursor_target_mz: Float64Builder::new(),
            precursor_selected_mz: Float64Builder::new(),
            precursor_isolation_width: Float64Builder::new(),
            precursor_charge: Int32Builder::new(),
            precursor_intensity: Float64Builder::new(),
            precursor_collision_energy: Float64Builder::new(),
            precursor_collision_energy_unit: StringBuilder::new(),
            precursor_native_id: StringBuilder::new(),
            precursor_activation: StringBuilder::new(),
            precursor_analyzer: StringBuilder::new(),
            precursor_ccs: Float64Builder::new(),
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
            acquisition_event_id: UInt32Builder::new(),
            // Field names are pinned so the builder always matches the
            // schema, whatever arrow's `MapBuilder` defaults are.
            extra: MapBuilder::new(
                Some(MapFieldNames {
                    entry: EXTRA_ENTRIES_FIELD.to_string(),
                    key: EXTRA_KEY_FIELD.to_string(),
                    value: EXTRA_VALUE_FIELD.to_string(),
                }),
                StringBuilder::new(),
                StringBuilder::new(),
            ),
            precursor_collision_energies: struct_list_builder(
                energy_fields(),
                vec![
                    Box::new(Float64Builder::new()),
                    Box::new(StringBuilder::new()),
                ],
            ),
            scan_windows: struct_list_builder(
                window_fields(),
                vec![
                    Box::new(Float64Builder::new()),
                    Box::new(Float64Builder::new()),
                ],
            ),
            warnings: struct_list_builder(
                warning_fields(),
                vec![
                    Box::new(StringBuilder::new()),
                    Box::new(StringBuilder::new()),
                ],
            ),
            seen_ids: HashSet::new(),
        }
    }

    /// Schema for the batch produced by this builder.
    pub fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    /// Append one spectrum row after checking it, so a malformed record is
    /// an error instead of a row that disagrees with itself. Nothing is
    /// appended when this returns `Err`.
    ///
    /// # Errors
    ///
    /// [`BatchError::Record`] when [`SpectrumRecord::validate`] fails or the
    /// native id was already pushed to this builder (across batches taken
    /// with [`finish_batch`](Self::finish_batch)), and
    /// [`BatchError::ColumnRange`] when `index` does not fit `u32` or
    /// `ms_level` does not fit `u8`.
    pub fn push(&mut self, rec: &SpectrumRecord) -> Result<(), BatchError> {
        let record = |problem| BatchError::Record {
            record_id: rec.native_id.clone(),
            problem,
        };
        rec.validate().map_err(record)?;
        if self.seen_ids.contains(&rec.native_id) {
            return Err(record(RecordProblem::DuplicateId));
        }
        let index = u32::try_from(rec.index).map_err(|_| BatchError::ColumnRange {
            record_id: rec.native_id.clone(),
            column: "index",
            value: rec.index as u64,
        })?;
        let ms_level = u8::try_from(rec.ms_level).map_err(|_| BatchError::ColumnRange {
            record_id: rec.native_id.clone(),
            column: "ms_level",
            value: u64::from(rec.ms_level),
        })?;
        self.seen_ids.insert(rec.native_id.clone());
        self.index.append_value(index);
        self.scan_number.append_value(rec.scan_number);
        self.native_id.append_value(&rec.native_id);
        self.ms_level.append_value(ms_level);
        self.polarity
            .append_option(rec.polarity.map(Polarity::as_str));
        self.scan_mode
            .append_option(rec.scan_mode.map(ScanMode::as_str));
        self.analyzer
            .append_option(rec.analyzer.map(Analyzer::as_str));
        self.filter.append_option(rec.filter.as_deref());
        self.retention_time_sec
            .append_option(rec.retention_time_sec);
        self.total_ion_current.append_option(rec.total_ion_current);
        self.base_peak_mz.append_option(rec.base_peak_mz);
        self.base_peak_intensity
            .append_option(rec.base_peak_intensity);
        self.low_mz.append_option(rec.low_mz);
        self.high_mz.append_option(rec.high_mz);
        self.ion_injection_time_ms
            .append_option(rec.ion_injection_time_ms);
        self.inv_mobility.append_option(rec.inv_mobility);
        self.faims_cv.append_option(rec.faims_cv);

        match &rec.precursor {
            Some(p) => {
                self.precursor_target_mz.append_option(p.target_mz);
                self.precursor_selected_mz.append_option(p.selected_mz);
                self.precursor_isolation_width
                    .append_option(p.isolation_width);
                self.precursor_charge.append_option(p.charge);
                self.precursor_intensity.append_option(p.intensity);
                self.precursor_collision_energy
                    .append_option(p.collision_energies.first().map(|ce| ce.value));
                self.precursor_collision_energy_unit.append_option(
                    p.collision_energies
                        .first()
                        .and_then(|ce| ce.unit)
                        .map(CollisionEnergyUnit::as_str),
                );
                for ce in &p.collision_energies {
                    let sb = self.precursor_collision_energies.values();
                    sb.field_builder::<Float64Builder>(0)
                        .expect("value builder")
                        .append_value(ce.value);
                    sb.field_builder::<StringBuilder>(1)
                        .expect("unit builder")
                        .append_option(ce.unit.map(CollisionEnergyUnit::as_str));
                    sb.append(true);
                }
                self.precursor_collision_energies.append(true);
                self.precursor_native_id
                    .append_option(p.precursor_native_id.as_deref());
                self.precursor_activation
                    .append_option(p.activation.map(Activation::as_str));
                self.precursor_analyzer
                    .append_option(p.analyzer.map(Analyzer::as_str));
                self.precursor_ccs.append_option(p.ccs);
            }
            None => {
                self.precursor_target_mz.append_null();
                self.precursor_selected_mz.append_null();
                self.precursor_isolation_width.append_null();
                self.precursor_charge.append_null();
                self.precursor_intensity.append_null();
                self.precursor_collision_energy.append_null();
                self.precursor_collision_energy_unit.append_null();
                self.precursor_collision_energies.append(false);
                self.precursor_native_id.append_null();
                self.precursor_activation.append_null();
                self.precursor_analyzer.append_null();
                self.precursor_ccs.append_null();
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
            .append_option(self.mobility_kind.map(MobilityArrayKind::as_str));
        self.acquisition_event_id
            .append_option(rec.acquisition_event_id);
        for (key, value) in &rec.extra {
            self.extra.keys().append_value(key);
            self.extra.values().append_value(value);
        }
        self.extra
            .append(true)
            .expect("extra keys and values are appended together");
        for w in &rec.scan_windows {
            let sb = self.scan_windows.values();
            sb.field_builder::<Float64Builder>(0)
                .expect("lower builder")
                .append_value(w.lower);
            sb.field_builder::<Float64Builder>(1)
                .expect("upper builder")
                .append_value(w.upper);
            sb.append(true);
        }
        self.scan_windows.append(true);
        for w in &rec.warnings {
            let sb = self.warnings.values();
            sb.field_builder::<StringBuilder>(0)
                .expect("code builder")
                .append_value(w.kind.code());
            sb.field_builder::<StringBuilder>(1)
                .expect("message builder")
                .append_value(&w.message);
            sb.append(true);
        }
        self.warnings.append(true);
        Ok(())
    }

    /// Number of rows accumulated so far.
    pub fn len(&self) -> usize {
        self.index.len()
    }

    /// `true` if no rows have been pushed.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Materialize the accumulated rows into a [`RecordBatch`] and consume
    /// the builder.
    ///
    /// # Errors
    ///
    /// [`BatchError::Arrow`] if Arrow rejects the arrays.
    pub fn finish(mut self) -> Result<RecordBatch, BatchError> {
        self.finish_batch()
    }

    /// Materialize the rows accumulated so far into a [`RecordBatch`] and
    /// keep the builder, including the native ids seen so far, so the next
    /// batch of a stream is still checked for duplicates against this one.
    ///
    /// # Errors
    ///
    /// [`BatchError::Arrow`] if Arrow rejects the arrays.
    pub fn finish_batch(&mut self) -> Result<RecordBatch, BatchError> {
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
            Arc::new(self.faims_cv.finish()),
            Arc::new(self.precursor_target_mz.finish()),
            Arc::new(self.precursor_selected_mz.finish()),
            Arc::new(self.precursor_isolation_width.finish()),
            Arc::new(self.precursor_charge.finish()),
            Arc::new(self.precursor_intensity.finish()),
            Arc::new(self.precursor_collision_energy.finish()),
            Arc::new(self.precursor_collision_energy_unit.finish()),
            Arc::new(self.precursor_native_id.finish()),
            Arc::new(self.precursor_activation.finish()),
            Arc::new(self.precursor_analyzer.finish()),
            Arc::new(self.precursor_ccs.finish()),
            Arc::new(self.mz.finish()),
            Arc::new(self.intensity.finish()),
            Arc::new(self.inv_mobility_per_peak.finish()),
            Arc::new(self.mobility_array_kind_col.finish()),
            Arc::new(self.acquisition_event_id.finish()),
            Arc::new(self.extra.finish()),
            Arc::new(self.precursor_collision_energies.finish()),
            Arc::new(self.scan_windows.finish()),
            Arc::new(self.warnings.finish()),
        ];
        Ok(RecordBatch::try_new(self.schema.clone(), arrays)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow_array::Array;
    use openmassspec_core::PrecursorInfo;

    fn rec(index: usize, ms_level: u32, n_peaks: usize, with_mob: bool) -> SpectrumRecord {
        let mz: Vec<f64> = (0..n_peaks).map(|i| 100.0 + i as f64).collect();
        let intensity: Vec<f32> = (0..n_peaks).map(|i| 10.0 + i as f32).collect();
        let mobility = if with_mob {
            Some((0..n_peaks).map(|i| 0.5 + i as f32 * 0.01).collect())
        } else {
            None
        };
        let mut rec = SpectrumRecord::new(index, (index + 1) as u32, format!("scan={}", index + 1))
            .with_ms_level(ms_level)
            .with_polarity(Polarity::Positive)
            .with_scan_mode(ScanMode::Centroid)
            .with_analyzer(Analyzer::TOFMS)
            .with_acquisition_event_id((index == 0).then_some(7))
            .with_retention_time_sec(index as f64)
            .with_total_ion_current(intensity.iter().map(|&v| v as f64).sum::<f64>())
            .with_base_peak_mz(mz.last().copied())
            .with_base_peak_intensity(intensity.last().map(|&v| v as f64))
            .with_low_mz(mz.first().copied())
            .with_high_mz(mz.last().copied())
            .with_faims_cv(if index == 0 { Some(-40.0) } else { None })
            .with_inv_mobility_per_peak(mobility)
            .with_peaks(mz, intensity);
        if ms_level >= 2 {
            rec.precursor = Some(
                PrecursorInfo::default()
                    .with_target_mz(500.0)
                    .with_selected_mz(500.5)
                    .with_isolation_width(2.0)
                    .with_charge(2),
            );
        }
        if index == 0 {
            rec = rec.with_extra("opentfraw.resolution", "120000");
        }
        rec
    }

    #[test]
    fn schema_round_trip() {
        let mut b = SpectrumBatchBuilder::new(Some(MobilityArrayKind::DriftTimeMilliseconds));
        b.push(&rec(0, 1, 3, true)).unwrap();
        b.push(&rec(1, 2, 4, false)).unwrap();
        let batch = b.finish().unwrap();
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.schema().fields().len(), 37);
        assert_eq!(
            check_schema_version(&batch.schema()).unwrap(),
            SCHEMA_VERSION
        );
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
        let faims_col = batch
            .column_by_name("faims_cv")
            .unwrap()
            .as_any()
            .downcast_ref::<arrow_array::Float64Array>()
            .unwrap();
        assert_eq!(faims_col.value(0), -40.0);
        assert!(faims_col.is_null(1));
        let event_col = batch
            .column_by_name("acquisition_event_id")
            .unwrap()
            .as_any()
            .downcast_ref::<arrow_array::UInt32Array>()
            .unwrap();
        assert_eq!(event_col.value(0), 7);
        assert!(event_col.is_null(1));
        let extra_col = batch
            .column_by_name("extra")
            .unwrap()
            .as_any()
            .downcast_ref::<arrow_array::MapArray>()
            .unwrap();
        assert_eq!(extra_col.value_length(0), 1);
        assert_eq!(extra_col.value_length(1), 0);
        let keys = extra_col
            .keys()
            .as_any()
            .downcast_ref::<arrow_array::StringArray>()
            .unwrap();
        let values = extra_col
            .values()
            .as_any()
            .downcast_ref::<arrow_array::StringArray>()
            .unwrap();
        assert_eq!(keys.value(0), "opentfraw.resolution");
        assert_eq!(values.value(0), "120000");
    }

    fn unwrap_record(r: Result<(), BatchError>) -> RecordProblem {
        match r.unwrap_err() {
            BatchError::Record { problem, .. } => problem,
            other => panic!("expected BatchError::Record, got {other:?}"),
        }
    }

    #[test]
    fn push_rejects_malformed_records_with_typed_errors_and_appends_nothing() {
        let mut b = SpectrumBatchBuilder::new(None);
        let mut short = rec(0, 1, 3, false);
        short.intensity.pop();
        assert!(matches!(
            unwrap_record(b.push(&short)),
            RecordProblem::ArrayLengthMismatch { .. }
        ));
        let mut nan = rec(0, 1, 3, false);
        nan.intensity[1] = f32::NAN;
        assert!(matches!(
            unwrap_record(b.push(&nan)),
            RecordProblem::NonFinite { .. }
        ));
        let mut unsorted = rec(0, 1, 3, false);
        unsorted.mz.swap(0, 2);
        assert!(matches!(
            unwrap_record(b.push(&unsorted)),
            RecordProblem::UnsortedMz { .. }
        ));
        let mut mob = rec(0, 1, 3, true);
        mob.inv_mobility_per_peak.as_mut().unwrap().pop();
        assert!(matches!(
            unwrap_record(b.push(&mob)),
            RecordProblem::ArrayLengthMismatch { .. }
        ));
        assert!(matches!(
            b.push(&rec(0, 256, 3, false)),
            Err(BatchError::ColumnRange {
                column: "ms_level",
                ..
            })
        ));
        assert!(b.is_empty());
        b.push(&rec(0, 1, 3, false)).unwrap();
        assert_eq!(b.len(), 1);
    }

    #[test]
    fn duplicate_ids_are_rejected_across_batches() {
        let mut b = SpectrumBatchBuilder::new(None);
        b.push(&rec(0, 1, 2, false)).unwrap();
        assert_eq!(
            unwrap_record(b.push(&rec(0, 1, 2, false))),
            RecordProblem::DuplicateId
        );
        let first = b.finish_batch().unwrap();
        assert_eq!(first.num_rows(), 1);
        assert!(b.is_empty());
        // The next batch still knows scan=1 was used.
        assert_eq!(
            unwrap_record(b.push(&rec(0, 1, 2, false))),
            RecordProblem::DuplicateId
        );
        b.push(&rec(1, 1, 2, false)).unwrap();
        assert_eq!(b.finish().unwrap().num_rows(), 1);
    }

    #[test]
    fn schema_version_is_checked() {
        let schema = spectrum_record_schema();
        assert_eq!(check_schema_version(&schema).unwrap(), SCHEMA_VERSION);
        let future = Schema::new(Vec::<Field>::new()).with_metadata(HashMap::from([(
            SCHEMA_VERSION_KEY.to_string(),
            (SCHEMA_VERSION + 1).to_string(),
        )]));
        assert!(matches!(
            check_schema_version(&future),
            Err(BatchError::SchemaVersion { .. })
        ));
        assert!(check_schema_version(&Schema::new(Vec::<Field>::new())).is_err());
    }

    #[test]
    fn new_columns_carry_units_windows_and_warnings() {
        use arrow_array::{Array, LargeListArray, StringArray, StructArray};
        use openmassspec_core::{Warning, WarningKind};
        let mut r = rec(1, 2, 2, false)
            .with_scan_window(100.0, 200.0)
            .with_scan_window(300.0, 400.0)
            .with_warning(Warning::new(WarningKind::UnknownScanMode, "no mode"));
        r.precursor = Some(
            PrecursorInfo::default()
                .with_selected_mz(500.0)
                .with_collision_energy(20.0, CollisionEnergyUnit::Nce)
                .with_collision_energy(35.0, None),
        );
        let mut b = SpectrumBatchBuilder::new(None);
        b.push(&r).unwrap();
        b.push(&rec(0, 1, 2, false)).unwrap();
        let batch = b.finish().unwrap();

        let unit = batch
            .column_by_name("precursor_collision_energy_unit")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(unit.value(0), "nce");
        assert!(unit.is_null(1));
        let list = |name: &str| {
            batch
                .column_by_name(name)
                .unwrap()
                .as_any()
                .downcast_ref::<LargeListArray>()
                .unwrap()
                .clone()
        };
        let energies = list("precursor_collision_energies");
        assert_eq!(energies.value_length(0), 2);
        assert!(energies.is_null(1));
        let rows = energies.value(0);
        let rows = rows.as_any().downcast_ref::<StructArray>().unwrap();
        let units = rows
            .column(1)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(units.value(0), "nce");
        assert!(units.is_null(1));
        let windows = list("scan_windows");
        assert_eq!(windows.value_length(0), 2);
        assert_eq!(windows.value_length(1), 0);
        let warnings = list("warnings");
        let w = warnings.value(0);
        let w = w.as_any().downcast_ref::<StructArray>().unwrap();
        let code = w.column(0).as_any().downcast_ref::<StringArray>().unwrap();
        let message = w.column(1).as_any().downcast_ref::<StringArray>().unwrap();
        assert_eq!(
            (code.value(0), message.value(0)),
            ("unknown_scan_mode", "no mode")
        );
    }
}
