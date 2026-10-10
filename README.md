# openmassspec-core

[![CI](https://github.com/Sigilweaver/OpenMassSpecCore/actions/workflows/ci.yml/badge.svg)](https://github.com/Sigilweaver/OpenMassSpecCore/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/openmassspec-core.svg)](https://crates.io/crates/openmassspec-core)
[![docs.rs](https://img.shields.io/docsrs/openmassspec-core)](https://docs.rs/openmassspec-core)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Rust MSRV](https://img.shields.io/badge/rust-1.88%2B-orange.svg)](https://www.rust-lang.org)

> Part of the [OpenMassSpec](https://sigilweaver.app/openmassspec/docs/)
> stack for mass spectrometry raw-file access. Sibling readers:
> [OpenTFRaw](https://github.com/Sigilweaver/OpenTFRaw) (Thermo),
> [OpenTimsTDF](https://github.com/Sigilweaver/OpenTimsTDF) (Bruker),
> [OpenWRaw](https://github.com/Sigilweaver/OpenWRaw) (Waters),
> [OpenARaw](https://github.com/Sigilweaver/OpenARaw) (Agilent),
> [OpenSXRaw](https://github.com/Sigilweaver/OpenSXRaw) (SCIEX),
> [OpenSZRaw](https://github.com/Sigilweaver/OpenSZRaw) (Shimadzu).

Shared, vendor-neutral foundation for the OpenMassSpec stack:
the `SpectrumSource` trait every parser implements, the canonical
`SpectrumRecord` / `RunMetadata` types, a streaming mzML 1.1.0 writer
(with optional indexed mzML output and SHA-1 footer), typed PSI-MS terms
generated from `psi-ms.obo`, and a cross-vendor conformance harness.

This repository is a Cargo workspace with two crates:

| Crate                | Purpose                                                        |
| -------------------- | -------------------------------------------------------------- |
| `openmassspec-core`  | Types, `SpectrumSource`, mzML writer, CV terms, conformance.   |
| `openmassspec-arrow` | Apache Arrow `RecordBatch` bridge for `SpectrumRecord` streams. |

- MSRV: 1.88
- License: Apache-2.0
- `#![forbid(unsafe_code)]`

Documentation: [sigilweaver.app/openmassspec/docs](https://sigilweaver.app/openmassspec/docs)

Reader authors upgrading from 2.x: see [MIGRATING-3.0.md](MIGRATING-3.0.md).

## Install

```sh
cargo add openmassspec-core
```

With the Arrow bridge:

```sh
cargo add openmassspec-arrow
```

## Quick example

Implement `SpectrumSource` and write a valid indexed mzML document. Decode
failures are yielded as `Err` items; the writer stops at the first one
unless you opt into skipping.

```rust,no_run
use openmassspec_core::{
    write_indexed_mzml, CvTerm, RunMetadata, SpectrumIter, SpectrumRecord, SpectrumSource,
};

struct MySource {
    spectra: Vec<SpectrumRecord>,
}

impl SpectrumSource for MySource {
    fn run_metadata(&self) -> RunMetadata {
        RunMetadata::new(
            "run.raw",
            CvTerm::new("MS:1000563", "Thermo RAW format"),
            CvTerm::new("MS:1000768", "Thermo nativeID format"),
            CvTerm::new("MS:1000483", "Thermo Fisher Scientific instrument model"),
            "my-reader",
            "1.0.0",
        )
    }
    fn iter_spectra<'a>(&'a mut self) -> SpectrumIter<'a> {
        Box::new(self.spectra.drain(..).map(Ok))
    }
    fn spectrum_count_hint(&self) -> Option<usize> {
        Some(self.spectra.len())
    }
}

let spectrum = SpectrumRecord::new(0, 1, "controllerType=0 controllerNumber=1 scan=1")
    .with_retention_time_sec(1.5)
    .with_peaks(vec![100.0, 200.0], vec![10.0, 20.0]);
let mut src = MySource { spectra: vec![spectrum] };
let mut out = std::fs::File::create("run.mzML")?;
let summary = write_indexed_mzml(&mut src, &mut out)?;
assert_eq!(summary.spectra_written, 1);
# Ok::<(), Box<dyn std::error::Error>>(())
```

To keep writing past undecodable records, pass
`WriteOptions::new().with_record_error_policy(RecordErrorPolicy::Skip)` to
`write_mzml_with_options` / `write_indexed_mzml_with_options`. The returned
`WriteSummary` reports how many spectra and chromatograms were skipped.

## Public API

| Symbol                                                      | Purpose                                                                     |
| ----------------------------------------------------------- | --------------------------------------------------------------------------- |
| `SpectrumRecord`                                            | Decoded spectrum: id, ms level, polarity, rt, peaks, precursor.             |
| `PrecursorInfo`                                             | Selected / isolated precursor, charge, activation, collision energy.        |
| `ChromatogramRecord`                                        | TIC / BPC / SRM trace.                                                      |
| `RunMetadata`                                               | Run-level CV terms: instrument, source format, native id format.            |
| `CvTerm`                                                    | A controlled-vocabulary term supplied by a reader.                          |
| `cv`                                                        | Typed PSI-MS / UO terms generated from `psi-ms.obo` (`cv::PEAK_PICKING`).   |
| `Polarity`, `Analyzer`, `ScanMode`, `MsPower`, `Activation` | Standard enumerations.                                                      |
| `MobilityArrayKind`                                         | Per-peak inverse-mobility / drift-time array kind.                          |
| `SpectrumSource`                                            | Trait every parser implements. Iterators yield `Result` items.              |
| `write_mzml`, `write_mzml_with_options`                     | Stream a `SpectrumSource` to a plain mzML 1.1.0 document.                   |
| `write_indexed_mzml`, `write_indexed_mzml_with_options`     | Same, with `<indexList>` + SHA-1 footer for byte-offset indexing.           |
| `WriteOptions`, `RecordErrorPolicy`, `WriteSummary`         | Writer codec and error policy; per-write counts of written/skipped records. |
| `Centroided`                                                | Streaming profile-to-centroid adapter for any `SpectrumSource`.             |
| `conformance::assert_source_invariants`                     | Check a live `SpectrumSource` for cross-vendor invariants.                  |
| `conformance::assert_iter_invariants`                       | Same, from any `IntoIterator<Item = SpectrumRecord>`.                       |
| `Error`                                                     | Aggregate `thiserror`-based error type.                                     |
| `openmassspec_arrow::SpectrumBatchBuilder`                  | Builder for `arrow_array::RecordBatch` from a spectrum stream.              |
| `openmassspec_arrow::spectrum_record_schema`                | The canonical Arrow schema.                                                 |

Public structs and enums are `#[non_exhaustive]`: build records with
constructors and `with_*` setters, and match enums with a wildcard arm.
Fields and variants can then be added without a major release.

## CV terms

Every PSI-MS / UO term the writer emits on its own is a typed constant in
`openmassspec_core::cv`, generated by `scripts/gen-cv.py` from a pinned
HUPO-PSI `psi-ms.obo` release (the source URL and `data-version` are in the
header of `src/cv/terms.rs` and in the emitted `<cv id="MS">` entry). The
build never fetches the ontology. To add a term, add a line to
`scripts/cv-terms.txt`, rerun the generator, then rerun
`scripts/extract-cv-subset.sh` so `tests/cv_terms.rs` can check it.

## Conformance harness

The conformance module enforces the cross-vendor invariants every parser
must satisfy:

- the source yields no `Err` items,
- monotonic spectrum indices,
- non-decreasing retention times per acquisition stream,
- equal-length m/z and intensity arrays,
- equal-length mobility arrays (when present),
- declared TIC and base-peak intensity match the arrays,
- precursor presence on MSn spectra.

Failures surface as `ConformanceError` variants
(`Source`, `PeakArrayLengthMismatch`, `MobilityArrayLengthMismatch`,
`RetentionTimeNonMonotonic`, and others).

The `vendor2mzml validate` subcommand in the
[OpenMassSpec](https://github.com/Sigilweaver/OpenMassSpec) umbrella runs
this harness on any vendor input or pre-existing mzML.

## Changelog

See [CHANGELOG.md](CHANGELOG.md).

## License

Apache-2.0. See [LICENSE](LICENSE).
