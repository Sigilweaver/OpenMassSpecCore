# openproteo-core

[![CI](https://github.com/Sigilweaver/OpenProteoCore/actions/workflows/ci.yml/badge.svg)](https://github.com/Sigilweaver/OpenProteoCore/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/openproteo-core.svg)](https://crates.io/crates/openproteo-core)
[![docs.rs](https://img.shields.io/docsrs/openproteo-core)](https://docs.rs/openproteo-core)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Rust MSRV](https://img.shields.io/badge/rust-1.85%2B-orange.svg)](https://www.rust-lang.org)

> Part of the [OpenProteo](https://sigilweaver.app/openproteo/docs/)
> stack for proteomics raw-file access. Sibling readers:
> [OpenTFRaw](https://github.com/Sigilweaver/OpenTFRaw) (Thermo),
> [OpenTimsTDF](https://github.com/Sigilweaver/OpenTimsTDF) (Bruker),
> [OpenWRaw](https://github.com/Sigilweaver/OpenWRaw) (Waters).

Shared, vendor-neutral foundation for the OpenProteo mass-spec stack:
the `SpectrumSource` trait every parser implements, the canonical
`SpectrumRecord` / `RunMetadata` types, a streaming mzML 1.1.0 writer
(with optional indexed mzML output and SHA-1 footer), an optional
Apache Arrow `RecordBatch` bridge, and a cross-vendor conformance
harness.

- MSRV: 1.85
- License: Apache-2.0
- `#![forbid(unsafe_code)]`

Documentation: [sigilweaver.app/openproteo/docs](https://sigilweaver.app/openproteo/docs)

## Install

```sh
cargo add openproteo-core
```

With the optional Arrow bridge:

```sh
cargo add openproteo-core --features arrow
```

## Quick example

Implement `SpectrumSource` and write a valid indexed mzML document:

```rust
use openproteo_core::{
    write_indexed_mzml, RunMetadata, SpectrumRecord, SpectrumSource,
};

struct MySource {
    spectra: Vec<SpectrumRecord>,
}

impl SpectrumSource for MySource {
    fn run_metadata(&self) -> RunMetadata {
        RunMetadata::default()
    }
    fn iter_spectra<'a>(&'a mut self) -> Box<dyn Iterator<Item = SpectrumRecord> + 'a> {
        Box::new(self.spectra.drain(..))
    }
    fn spectrum_count(&self) -> Option<usize> {
        Some(self.spectra.len())
    }
}

let mut src = MySource { spectra: vec![/* SpectrumRecord { .. } */] };
let mut out = std::fs::File::create("run.mzML")?;
write_indexed_mzml(&mut src, &mut out)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Public API

| Symbol                                                      | Module          | Purpose                                                                  |
| ----------------------------------------------------------- | --------------- | ------------------------------------------------------------------------ |
| `SpectrumRecord`                                            | `types`         | Decoded spectrum: id, ms level, polarity, rt, peaks, precursor.          |
| `PrecursorInfo`                                             | `types`         | Selected / isolated precursor, charge, activation, scan window.          |
| `ChromatogramRecord`                                        | `types`         | TIC / BPC / SRM trace.                                                   |
| `RunMetadata`                                               | `types`         | Run-level CV terms: instrument, source format, native id format.         |
| `CvTerm`                                                    | `types`         | A PSI-MS controlled-vocabulary term.                                     |
| `Polarity`, `Analyzer`, `ScanMode`, `MsPower`, `Activation` | `enums`         | Standard enumerations.                                                   |
| `MobilityArrayKind`                                         | `enums`         | Per-peak inverse-mobility / drift-time array kind.                       |
| `SpectrumSource`                                            | `source`        | Trait every parser implements.                                           |
| `write_mzml`                                                | `mzml`          | Stream a `SpectrumSource` to a plain mzML 1.1.0 document.                |
| `write_indexed_mzml`                                        | `mzml`          | Same, with `<indexList>` + SHA-1 footer for byte-offset indexing.        |
| `conformance::assert_source_invariants`                     | `conformance`   | Check a live `SpectrumSource` for cross-vendor invariants.               |
| `conformance::assert_iter_invariants`                       | `conformance`   | Same, but from any `IntoIterator<Item = SpectrumRecord>`.                |
| `arrow::SpectrumBatchBuilder`                               | `arrow` (feat)  | Zero-copy builder for `arrow_array::RecordBatch` from a spectrum stream. |
| `arrow::spectrum_record_schema`                             | `arrow` (feat)  | The canonical Arrow schema.                                              |
| `Error`                                                     | `error`         | Aggregate `thiserror`-based error type.                                  |

## Conformance harness

The conformance module enforces the cross-vendor invariants every parser
must satisfy:

- monotonic spectrum indices,
- non-negative, non-decreasing retention times,
- equal-length m/z and intensity arrays,
- equal-length mobility arrays (when present),
- MS-level / polarity sanity,
- precursor presence on MSn spectra.

Failures surface as `ConformanceError` variants
(`PeakArrayLengthMismatch`, `MobilityArrayLengthMismatch`,
`RetentionTimeNonMonotonic`, and others).

The `vendor2mzml validate` subcommand in the
[OpenProteo](https://github.com/Sigilweaver/OpenProteo) umbrella runs
this harness on any vendor input or pre-existing mzML.

## Feature flags

| Flag    | Default | Effect                                                    |
| ------- | :-----: | --------------------------------------------------------- |
| `arrow` |   off   | Enables `arrow_array::RecordBatch` building from spectra. |

## Changelog

See [CHANGELOG.md](CHANGELOG.md).

## License

Apache-2.0. See [LICENSE](LICENSE).
