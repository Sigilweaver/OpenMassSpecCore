# openmassspec-arrow

Apache Arrow bridge for
[`openmassspec-core`](https://crates.io/crates/openmassspec-core) spectrum
records: a canonical Arrow schema (`spectrum_record_schema`) and a streaming
`SpectrumBatchBuilder` that turns `SpectrumRecord`s into one
`arrow_array::RecordBatch`, one row per spectrum.

```rust
use openmassspec_arrow::SpectrumBatchBuilder;
use openmassspec_core::SpectrumRecord;

let mut builder = SpectrumBatchBuilder::new(None);
builder
    .push(&SpectrumRecord::new(0, 1, "scan=1").with_peaks(vec![100.0], vec![5.0]))
    .unwrap();
let batch = builder.finish().unwrap();
assert_eq!(batch.num_rows(), 1);
```

`push` validates each record with `SpectrumRecord::validate` and rejects a
repeated native id; the error is a `BatchError` that carries Core's
`RecordProblem`. The schema metadata holds `openmassspec.schema_version`;
`check_schema_version` rejects a newer one.

The arrow-rs release is chosen by a version-named feature: `arrow-60`
(default) uses arrow 60. This crate has its own version line. MSRV 1.88.
Apache-2.0.
