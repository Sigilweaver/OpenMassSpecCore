# Changelog

All notable changes to `openmassspec-core` are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
crate adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [1.4.0] - 2026-07-28

### Added

- The mzML writer now zlib-compresses (`MS:1000574`) binary data arrays
  (m/z, intensity, per-peak ion mobility, and chromatogram time/intensity)
  by default in `write_mzml`/`write_indexed_mzml`, instead of always
  emitting them uncompressed under a hardcoded `MS:1000576` ("no
  compression") cvParam. Output from downstream vendor crates was
  multiples larger than msconvert/OpenMS output for the same run as a
  result. New `write_mzml_with_compression`/
  `write_indexed_mzml_with_compression` entry points and a `Compression`
  enum (`Zlib` default, `NoCompression` to keep the previous byte layout)
  give callers explicit control; the `<fileChecksum>` SHA-1 in indexed
  output needed no code change since it already hashes the bytes actually
  written, whichever codec produced them. Numpress support is left as a
  follow-up. Contributed by @Nabejo (closes #8).
- Unit tests for the hand-rolled SHA-1 implementation backing indexed
  mzML's `<fileChecksum>` (`src/mzml.rs`), checked against NIST/FIPS
  180-1 test vectors and message lengths chosen to exercise the
  padding-block boundaries (empty, 56, 64, and 1000 bytes), plus a
  streaming-vs-single-shot equivalence check. Also adds an integration
  test that generates a small indexed-mzML document, recomputes its
  checksum with a second, independent SHA-1 implementation, and asserts
  it matches the embedded `fileChecksum`. Previously nothing in the
  crate would have caught a future refactor silently breaking the hash
  (contributed by @Nabejo, closes #11).
- `RunMetadata::acquisition_software_name` and
  `RunMetadata::acquisition_software_version`, both `Option<String>`,
  carrying the name/version of the original vendor acquisition software
  that recorded a run (Xcalibur, MassHunter, Analyst, LabSolutions, etc).
  This is distinct from the existing `software_name`/`software_version`
  fields, which identify the vendor-reader crate producing the mzML, not
  the instrument-control software that acquired the data. `softwareList`
  previously hardcoded a single `<software>` entry for the reader itself,
  so this provenance was dropped even where a vendor reader's own
  instrument-method metadata already exposes it - msconvert-produced
  files preserve it, ours did not. The writer now emits a second
  `<software>` entry (fixed id `acquisition_software`, since an
  arbitrary vendor-reported name may contain characters, such as a
  space, that are not valid in an XML `xs:ID`) when
  `acquisition_software_name` is present, using the same generic
  `MS:1000799` ("custom unreleased software tool") cvParam already used
  for the reader's own entry; mapping specific vendor names to their own
  PSI-MS CV term is left for a follow-up; each vendor string would need
  to be checked against the CV rather than guessed. Adding these fields
  is a source-breaking change for any code that constructs a
  `RunMetadata` literal - existing vendor crates need two lines added
  (both may default to `None`) at their `RunMetadata` construction
  site(s) (closes #10, contributed by @Nabejo).

### Fixed

- `write_mzml`/`write_indexed_mzml` now emit an `MS:1000035` ("peak
  picking") `<processingMethod>` step in the `<dataProcessingList>` when
  the spectrum stream passed in has been wrapped in `Centroided`;
  previously only the blanket `MS:1000544` ("Conversion to mzML") step
  was ever recorded, so mzML written from a centroided stream carried
  per-spectrum centroid-mode cvParams with no matching processing-history
  entry - a provenance gap most mzML consumers/validators expect to be
  internally consistent. Added `SpectrumSource::additional_processing_steps`
  (default: none) so adapters can advertise steps like this; sources that
  don't use `Centroided` are unaffected and produce byte-for-byte
  identical output (closes #9, contributed by @Nabejo).

## [1.3.0] - 2026-07-24

### Added

- `RunMetadata::instrument_serial_number`, an `Option<String>` carrying the
  physical instrument's serial number. The writer emits it as an
  `MS:1000529` ("instrument serial number") cvParam on the default
  `instrumentConfiguration` when present. Adding this field is a
  source-breaking change for any code that constructs a `RunMetadata`
  literal - existing vendor crates need one line added at their
  `RunMetadata` construction site(s) (closes #7).
- `RunMetadata::analyzers`, a `Vec<Analyzer>` declaring the instrument's
  physical mass-analyzer components. Each entry gets its own
  `<instrumentConfiguration>` (`IC2`, `IC3`, ...) with the matching PSI-MS
  analyzer-component cvParam, and `write_spectrum` now reads
  `SpectrumRecord::analyzer` (previously captured but never emitted) to add
  a per-spectrum `instrumentConfigurationRef`, so hybrid instruments that
  switch analyzers scan-to-scan (e.g. quad/ion-trap MS2 vs Orbitrap MS1 on
  a Fusion) keep that identity instead of losing it to a single blanket
  `IC1`. Leaving `analyzers` empty keeps the previous single-`IC1`-for-
  everything output byte-for-byte. Also source-breaking the same way as
  `instrument_serial_number` above (closes #6).
- `PrecursorInfo::ccs`, an `Option<f64>` carrying the selected ion's
  collision cross-sectional area in square angstroms. The writer emits it
  as an `MS:1002954` ("collisional cross sectional area") cvParam on the
  precursor's `selectedIon` when present, and it round-trips through the
  optional `arrow` feature's `precursor_ccs` column. This is the shared-
  schema field OpenTimsTDF#14 (Bruker 1/K0) and OpenWRaw#10 (Waters drift
  time) need once they convert their vendor-native ion-mobility values to
  CCS - previously there was nowhere in this crate to put the result.
  Source-breaking the same way as the two entries above (closes #5).

### Fixed

- The `arrow` feature's `SpectrumBatchBuilder` now includes a `faims_cv`
  column. `SpectrumRecord::faims_cv` (added in 1.2.0) was the one field
  missing from the schema's otherwise-complete 1:1 mirror of
  `SpectrumRecord`/`PrecursorInfo` - Arrow consumers silently lost FAIMS
  compensation voltage even though the mzML writer emitted it correctly.

## [1.2.0] - 2026-07-15

### Added

- `SpectrumRecord::faims_cv`, an `Option<f64>` carrying FAIMS compensation
  voltage in volts. The writer emits it as a scan-level `MS:1001581`
  ("FAIMS compensation voltage") cvParam when present. Adding this field
  is a source-breaking change for any code that constructs a
  `SpectrumRecord` literal without `..Default::default()` (there is no
  `Default` impl) - existing vendor crates need one line added at their
  `SpectrumRecord` construction site(s) (closes #3).

### Fixed

- `write_mzml`/`write_indexed_mzml` now call `SpectrumSource::iter_chromatograms`
  and emit a `<chromatogramList>` (with a second `<index name="chromatogram">`
  block in the indexed variant) when the source yields anything; previously
  chromatogram data had no path to output regardless of vendor (closes #1).
- `write_prologue` now emits the `<run startTimeStamp="...">` attribute
  from `RunMetadata.start_timestamp` when present. All five vendor crates
  already populate this field; the writer was silently dropping it
  (closes #2).

## [1.1.1] - 2026-07-12

### Changed

- Bumped `arrow-array`/`arrow-schema`/`arrow-buffer` from `^58` to `^59`
  (optional `arrow` feature). No public API change; this unblocks
  downstream crates that need `arrow`'s `pyarrow` feature to build
  against pyo3 0.29 (only one pyo3 version can be linked per binary).

## [1.1.0] - 2026-07-12

### Added

- `Centroided<S>`, a `SpectrumSource` adapter that centroids every
  profile-mode spectrum a wrapped source yields (local-maxima peak
  picking; the picked m/z, and inverse mobility when present, is the
  intensity-weighted centroid over each apex and its immediate
  neighbors). Spectra already tagged `ScanMode::Centroid` pass through
  unchanged, so wrapping a source is idempotent. An optional
  `with_min_intensity` builder method discards picked peaks below a
  noise floor. Composes with `write_mzml`/`write_indexed_mzml` and the
  Arrow bridge with no special-casing, since it is just another
  `SpectrumSource`.

## [1.0.0] - 2026-07-10

Renamed from `openproteo-core`. The vendor raw-file readers this crate
underpins (Thermo, Bruker, Waters, with Agilent and SCIEX joining the
suite) are used as much in metabolomics and lipidomics as in proteomics,
so the umbrella naming moved from proteomics-specific to general mass
spectrometry. No API or behavioral changes from `openproteo-core` 1.0.1;
version reset to 1.0.0 to reflect that this is a new package identity on
crates.io (the old `openproteo-core` name stays published and frozen at
1.0.1, it is not superseded in place). See
[openproteo-core's CHANGELOG](https://github.com/Sigilweaver/OpenProteoCore/blob/main/CHANGELOG.md)
for pre-rename history.

### Changed

- Package renamed `openproteo-core` -> `openmassspec-core`.

## [1.0.1] - 2026-05-22

Documentation polish to bring the crate landing page in line with the
rest of the OpenProteo stack. No API or behavioural changes.

### Changed

- README rewritten with CI / crates.io / docs.rs / license badges, a
  stack callout pointing at the sibling vendor readers, and a link to
  the unified docs hub at `sigilweaver.app/openproteo/docs`.
- `Cargo.toml`: `homepage` now points to the docs site, `documentation`
  field added (docs.rs), and a `[package.metadata.docs.rs]` block was
  added so docs.rs renders all features (`arrow`).

### Removed

- `ROADMAP.md` (internal planning artifact; no longer tracked).

## [1.0.0] - 2026-05-21

First stable release. No API changes from `0.1.0`; promoted to `1.0.0`
to align with the rest of the OpenProteo stack and to make the crate's
stability contract explicit. `0.1.0` has been yanked from crates.io.

### Changed

- MSRV bumped from 1.75 to 1.85 to track the `arrow-58.x` toolchain
  requirement (`edition = "2024"` Cargo feature) and to align with the
  rest of the OpenProteo stack.

## [0.1.0]

Initial published shape of the crate. This release defines the
vendor-neutral foundation the vendor parsers
(`opentfraw`, `opentimstdf`, `openwraw`) build on.

### Added

- Vendor-neutral record types: `SpectrumRecord`, `PrecursorInfo`,
  `ChromatogramRecord`, `RunMetadata`, `CvTerm`.
- Standard enumerations: `Polarity`, `Analyzer`, `ScanMode`, `MsPower`,
  `Activation`, `MobilityArrayKind`.
- `SpectrumSource` trait that every vendor parser implements; default
  empty `iter_chromatograms` and `spectrum_count`.
- Canonical mzML 1.1.0 writer (`write_mzml`) and indexed-mzML writer
  (`write_indexed_mzml`) with `<indexList>` and SHA-1 footer.
- Conformance harness (`assert_source_invariants` /
  `assert_iter_invariants`) with structured `ConformanceError`
  variants (peak-array length, mobility-array length, retention-time
  monotonicity, MS-level / polarity, precursor presence, index
  sequence, empty spectrum).
- Optional `arrow` feature: zero-copy `SpectrumBatchBuilder` and the
  canonical `spectrum_record_schema()` for downstream Arrow / Parquet
  / Lance consumers.
- Aggregate `Error` enum (`thiserror`-based) covering I/O, decode, and
  conformance failures.

### Policy

- MSRV pinned at Rust 1.75.
- `#![forbid(unsafe_code)]` crate-wide.
- License: Apache-2.0.

[1.0.0]: https://github.com/Sigilweaver/OpenMassSpecCore/releases/tag/v1.0.0-openmassspec
[0.1.0]: https://github.com/Sigilweaver/OpenMassSpecCore/releases/tag/v0.1.0
