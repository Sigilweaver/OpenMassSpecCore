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
builder.push(&SpectrumRecord::new(0, 1, "scan=1").with_peaks(vec![100.0], vec![5.0]));
let batch = builder.finish().unwrap();
assert_eq!(batch.num_rows(), 1);
```

Uses arrow 60. MSRV 1.88. Apache-2.0.
