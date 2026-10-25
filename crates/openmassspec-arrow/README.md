# openmassspec-arrow

Apache Arrow bridge for
[`openmassspec-core`](https://crates.io/crates/openmassspec-core) spectrum
records: a canonical Arrow schema (`spectrum_record_schema`) and a streaming
`SpectrumBatchBuilder` that turns `SpectrumRecord`s into one
`arrow_array::RecordBatch`, one row per spectrum.

```rust
use openmassspec_arrow::arrow_60::SpectrumBatchBuilder;
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

Each arrow-rs release has its own module and feature: `arrow_60` (feature
`arrow-60`, default) is built against arrow 60. A new arrow-rs major is a new
module beside it, generated from the same source, so existing paths keep
working; no arrow type appears at the crate root. This crate has its own
version line.

The schema metadata carries `openmassspec.schema_version`. Appending a column
does not change it; any other change to the schema, including a new field
inside a list-of-struct column, bumps it. A reader rejects a newer version and
a schema with no version, which includes files written by the 3.x crate.
MSRV 1.88. Apache-2.0.
