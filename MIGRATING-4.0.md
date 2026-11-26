# Migrating to openmassspec-core 4.0

What changes between `openmassspec-core` 3.0.0 and 4.0.0, and between
`openmassspec-arrow` 3.0.0 and 4.0.0, for three kinds of reader:

| You are                                                          | Read                                   |
| ---------------------------------------------------------------- | -------------------------------------- |
| Writing a Rust crate that implements `SpectrumSource`            | [1. Reader authors](#1-reader-authors) |
| Using `SpectrumRecord`, `PrecursorInfo` and the writer from Rust | [2. Users of records](#2-users-of-records) |
| Calling a reader from Python                                     | [3. Python users](#3-python-users)     |

The complete list of changes is the 4.0.0 entry of
[CHANGELOG.md](CHANGELOG.md). This note shows the edits each change asks of
you. Every `rust` block that follows an "After" heading is compiled as a doc
test by the unpublished `migrating-doctests` crate in this workspace, which
`cargo test --workspace` runs in CI; blocks marked "Before" show 3.0.0 code
and are not compiled.

Three changes decide most of the work:

- Collision energy is a list, and each entry names its unit or says it is
  unknown. Core no longer assumes eV.
- Retention time is optional, and a spectrum with no scan mode is an error.
  A record that Core 3 wrote with a silent default is now rejected.
- Iterators are `Send`, so a reader generic over `Read + Seek` adds `Send`
  to its bounds. That is a public signature change in the reader, so the
  reader needs a major release too.
- The mzML writer and `Centroided` check each record, and some output
  changes with no compile error. Section 1.5, 1.6 and 1.10 list them for
  reader authors and 2.3 for everyone who compares output.

## 1. Reader authors

Work through these in order. `cargo check` after each step points at what
is left.

### 1.1 Bump the dependency

```toml
[dependencies]
openmassspec-core = "4.0.0"
# Only if you use the Arrow bridge (see 1.9):
openmassspec-arrow = "4.0.0"
```

### 1.2 Iterators are `Send`

`SpectrumIter` and `ChromatogramIter` are `Box<dyn Iterator<..> + Send + 'a>`.
The iterator borrows `&mut self` and moves what it captures, so whatever it
captures must be `Send`. A reader generic over its byte source captures that
source, so it adds `Send` to the bound, on the struct and on the impl. A
captured value that is not `Send` (`Rc`, `RefCell`, a raw pointer) is
replaced by `Arc`, `Mutex` or a plain owned value, or kept out of the
iterator.

Before:

```rust,ignore
pub struct MyReader<R: Read + Seek> {
    file: R,
    scan_count: u32,
}

impl<R: Read + Seek> SpectrumSource for MyReader<R> {
    // ...
}
```

Without the change, rustc reports that `R` "cannot be sent between threads
safely" at the `Box::new` of the iterator and suggests restricting `R: Send`.

After:

```rust
use std::io::{Read, Seek, SeekFrom};

use openmassspec_core::{CvTerm, Result, RunMetadata, SpectrumIter, SpectrumRecord, SpectrumSource};

pub struct MyReader<R: Read + Seek + Send> {
    file: R,
    scan_count: u32,
}

impl<R: Read + Seek + Send> SpectrumSource for MyReader<R> {
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
        let file = &mut self.file;
        Box::new((0..self.scan_count).map(move |i| -> Result<SpectrumRecord> {
            let mut buf = [0u8; 16];
            file.seek(SeekFrom::Start(u64::from(i) * 16))?;
            file.read_exact(&mut buf)?;
            Ok(SpectrumRecord::new(i as usize, i + 1, format!("scan={}", i + 1)))
        }))
    }
}
```

If your public API exposes `R: Read + Seek`, the added `Send` makes your
crate's own release a major one.

### 1.3 Collision energy is a list with units

`PrecursorInfo::collision_energy` and `ce_is_nce` are gone, and so is
`with_ce_is_nce`. `PrecursorInfo::collision_energies` holds one
`CollisionEnergy { value, unit }` per energy, and
`with_collision_energy(value, unit)` appends one. There is no default unit.
Pass the unit the file states, or `None` when the file gives a value and
does not say what it is. Stepped or combined activation is several calls. A
`None` value appends nothing.

Before:

```rust,ignore
let precursor = PrecursorInfo::default()
    .with_collision_energy(30.0)
    .with_ce_is_nce(true);
```

After:

```rust
use openmassspec_core::{CollisionEnergyUnit, PrecursorInfo};

// The file says the energy is normalized (NCE) or in electronvolts.
let nce = PrecursorInfo::default().with_collision_energy(30.0, CollisionEnergyUnit::Nce);
let ev = PrecursorInfo::default().with_collision_energy(35.0, CollisionEnergyUnit::Ev);

// Stepped collision energy: one entry per step.
let stepped = PrecursorInfo::default()
    .with_collision_energy(20.0, CollisionEnergyUnit::Nce)
    .with_collision_energy(30.0, CollisionEnergyUnit::Nce);
assert_eq!(stepped.collision_energies.len(), 2);

// A value whose unit the file does not state: do not guess.
let unknown = PrecursorInfo::default().with_collision_energy(30.0, None::<CollisionEnergyUnit>);
assert_eq!(unknown.collision_energies[0].unit, None);

// No energy in the file: nothing is appended.
let none = PrecursorInfo::default().with_collision_energy(None::<f64>, CollisionEnergyUnit::Nce);
assert!(none.collision_energies.is_empty());
```

A reader whose 3.x release wrote `ce_is_nce: false` because it did not know
the unit now passes `None`. The mzML writer rejects a collision energy with
no unit (`RecordProblem::CollisionEnergyUnitUnknown`) unless the caller
allows `Check::UnknownCollisionEnergyUnit` (see 1.6), which writes the value
as a unit-less `userParam` and adds a warning to the record. A reader whose
format never states the unit declares that once for the whole source (see
1.6, "A format that never states a value").

An entry is one applied energy. A ramp is not an entry: do not put the
start, end or midpoint of a ramp in `collision_energies`, and do not put one
cell's value there when the scan uses two cells, because each states a
value the instrument did not apply. A reader that cannot state the applied
energy leaves the list empty and keeps the raw values in `extra`. Stepped
energies are different and fit: one entry per step. `CollisionEnergy` has no
ramp or cell-role field in 4.0; that is a known gap, and fields for it can
be added in a 4.x release.

### 1.4 Retention time is optional

`SpectrumRecord::retention_time_sec` is `Option<f64>`, and
`SpectrumRecord::new` leaves it `None` (it used to set `0.0`). The setter
takes seconds as an `f64` or an `Option<f64>` (`impl Into<Option<f64>>`).

A spectrum with no retention time is an error
(`RecordProblem::MissingRetentionTime`) unless the run declares that the
format has no time axis. In that case no spectrum may carry one
(`RecordProblem::UnexpectedRetentionTime`). A reader that never set the
field because `0.0` was the default must now set it from the file.

Before:

```rust,ignore
let rec = SpectrumRecord::new(0, 1, "scan=1")
    .with_retention_time_sec(meta.retention_time_min * 60.0);
```

After:

```rust
use openmassspec_core::{CvTerm, RunMetadata, SpectrumRecord, TimeAxis};

// The decoder may or may not have found a time: pass the Option through.
let decoded_minutes: Option<f64> = Some(12.5);
let rec = SpectrumRecord::new(0, 1, "scan=1")
    .with_retention_time_sec(decoded_minutes.map(|minutes| minutes * 60.0));
assert_eq!(rec.retention_time_sec, Some(750.0));

// A format with no time axis says so once, on the run.
let run = RunMetadata::new(
    "run.dat",
    CvTerm::new("MS:1000563", "Thermo RAW format"),
    CvTerm::new("MS:1000768", "Thermo nativeID format"),
    CvTerm::new("MS:1000483", "Thermo Fisher Scientific instrument model"),
    "my-reader",
    "1.0.0",
)
.with_time_axis(TimeAxis::Absent);
assert_eq!(run.time_axis, TimeAxis::Absent);
```

Code that reads the field handles `None` explicitly. `unwrap_or(0.0)` puts
back the silent default this change removes.

### 1.5 Scan mode and scan windows

Set `scan_mode` on every spectrum from what the file declares. A spectrum
with `scan_mode: None` is `RecordProblem::UnknownScanMode` and stops the
write even under `RecordErrorPolicy::Skip`. Core 3 wrote it as a profile
spectrum. A reader that cannot tell leaves `None` and lets the caller opt in
with `Check::UnknownScanMode` (1.6), or, when its format never states the
mode in any file, declares that in `RunMetadata::unstated` (1.6).

`Centroided` picks peaks only from spectra declared `ScanMode::Profile`, and
runs the record checks below first. It rejects a spectrum labelled profile
whose median m/z spacing is more than 0.5 percent of m/z
(`RecordProblem::NotProfile`), which is how a peak list mislabelled as
profile shows up: fix the label in the reader (`ScanMode::Centroid`) instead
of working around it. It joins neighbouring points into one peak only
within three median spacings, so a gap in zero-stripped profile data splits
peaks where Core 3 averaged across it.

`SpectrumRecord::scan_windows` is new: the instrument's acquisition ranges
in m/z, one entry per window. The writer emits `<scanWindowList>` only from
it. Core 3 wrote one window built from the observed m/z range, which is not
the instrument setting, so a reader that knows the acquisition range sets it
and one that does not writes no window.

```rust
use openmassspec_core::{Polarity, ScanMode, SpectrumRecord};

let rec = SpectrumRecord::new(0, 1, "scan=1")
    .with_polarity(Polarity::Positive)
    .with_scan_mode(ScanMode::Centroid)
    // One call per window; `with_scan_windows` takes a whole `Vec<ScanWindow>`.
    .with_scan_window(100.0, 1500.0)
    .with_peaks(vec![150.1, 250.2], vec![10.0, 20.0]);
assert_eq!(rec.scan_windows.len(), 1);
assert!(rec.validate().is_ok());
```

`SpectrumRecord::validate` runs the record checks that the writer, the Arrow
bridge and `Centroided` share. Besides the ones above and those in 1.3 and
1.4, a reader will meet these:

| Trigger                                                                                           | Error                                              | What the reader does                                                                                       |
| ------------------------------------------------------------------------------------------------- | -------------------------------------------------- | ---------------------------------------------------------------------------------------------------------- |
| A scan window with a negative lower limit, `lower > upper`, or a non-finite limit                  | `RecordProblem::InvalidScanWindow` (`NonFinite`)   | Fix the decoder; add no window when the file states none                                                  |
| Precursor `charge`, `intensity` or `ccs` set while `selected_mz` is `None`                         | `RecordProblem::PrecursorWithoutSelectedMz`        | Set `selected_mz` from the file, or leave the dependent field unset                                       |
| Spectrum-level `inv_mobility` set while the run's `mobility_array_kind` is `DriftTimeMilliseconds` | `RecordProblem::ScanMobilityKindUnsupported`       | Leave `inv_mobility` unset; carry drift times per peak in `inv_mobility_per_peak`                         |
| `intensity` or `inv_mobility_per_peak` not the length of `mz`                                      | `RecordProblem::ArrayLengthMismatch`               | Fix the decoder; the mobility array used to be dropped silently                                           |
| NaN or infinite values, an `mz` array that is not ascending, `ms_level` 0                          | `RecordProblem::NonFinite`, `UnsortedMz`, `MsLevelZero` | Emit sorted arrays from the decoder (see below)                                                        |
| A chromatogram whose `intensity` is not the length of `time_sec`                                   | `RecordProblem::ArrayLengthMismatch`               | Fix the decoder                                                                                            |
| NaN or infinite values in a chromatogram's `time_sec`, `intensity`, `precursor_mz` or `product_mz` | `RecordProblem::NonFinite`                         | Reject or repair in the decoder                                                                            |
| Two chromatograms with the same `id`                                                               | `RecordProblem::DuplicateId`                       | Give each chromatogram a distinct id                                                                       |

`ChromatogramRecord::validate` runs the chromatogram length and finite-value
checks; the writer adds the repeated-id check. No `Check` covers any of
them. In particular unsorted m/z cannot be allowed: no
`Leniency` and no `--allow` name relaxes it, so a reader must emit `mz`
sorted ascending, with `intensity` and `inv_mobility_per_peak` reordered
with it, and must not rely on the caller to opt in. The same holds for
non-finite values, mismatched array lengths and `ms_level` 0.

The scan-level drift-time row affects readers of drift-time mobility data
(Waters-style): `RunMetadata::mobility_array_kind` is the run-wide kind, and
`None` means 1/K0.

### 1.6 Leniency, `Check` and warnings

By default the writer and `Centroided` reject what they cannot represent
faithfully. A `Leniency` names the `Check`s that may continue instead, and
each affected record gets a `Warning`:

| Check                          | Strict mode                                | Lenient mode                                    |
| ------------------------------ | ------------------------------------------ | ----------------------------------------------- |
| `UnknownScanMode`              | record error                               | written with no profile or centroid term        |
| `DanglingPrecursorRef`         | record error: the referenced spectrum is not among those already written | the `spectrumRef` is dropped          |
| `UnknownCollisionEnergyUnit`   | record error                               | value written as a unit-less `userParam`        |
| `SpectrumCountShortfall`       | `Error::SpectrumCountShortfall`            | run-level warning and `WriteSummary::spectrum_count_mismatch` |

Each error message names the check that allows it (for example
"allow the check `unknown_scan_mode` to continue without it"), so a
command line that exposes `--allow` can show the user what to pass.

Precursor references depend on stream order. The writer resolves
`PrecursorInfo::precursor_native_id` only against spectra it has already
written, so a precursor that names a spectrum later in the stream, or one
that was skipped or never yielded, is dangling. An MS2 spectrum yielded
before its MS1 is a `RecordProblem::DanglingPrecursorRef` even though the MS1
appears later in the output. A reader yields each MS1 before the spectra
that reference it, or leaves `precursor_native_id` unset when the parent is
not part of the stream.

```rust
use openmassspec_core::{
    write_mzml, CvTerm, Error, PrecursorInfo, RecordProblem, RunMetadata, ScanMode,
    SpectrumIter, SpectrumRecord, SpectrumSource,
};

struct Fixed(Vec<SpectrumRecord>);

impl SpectrumSource for Fixed {
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
        Box::new(self.0.clone().into_iter().map(Ok))
    }
}

fn spectrum(index: usize, id: &str) -> SpectrumRecord {
    SpectrumRecord::new(index, index as u32 + 1, id)
        .with_scan_mode(ScanMode::Centroid)
        .with_retention_time_sec(index as f64)
        .with_peaks(vec![100.0], vec![1.0])
}

let ms1 = spectrum(0, "scan=1");
let ms2 = spectrum(1, "scan=2")
    .with_ms_level(2)
    .with_precursor(PrecursorInfo::default().with_precursor_native_id("scan=1".to_string()));

// The MS2 comes first: its parent is not written yet.
let mut wrong_order = Fixed(vec![ms2.clone(), ms1.clone()]);
let err = write_mzml(&mut wrong_order, &mut Vec::new()).unwrap_err();
assert!(matches!(
    err,
    Error::Record { problem: RecordProblem::DanglingPrecursorRef { .. }, .. }
));

// Parent first: fine.
let mut right_order = Fixed(vec![ms1, ms2]);
write_mzml(&mut right_order, &mut Vec::new()).unwrap();
```

`Check` is open. A reader defines its own with `Check::other` and a
namespaced code, parses its own names, and attaches a warning with
`WarningKind::Other` to each record it let through. Core never consults a
reader-defined check: the reader receives the `Leniency` from its caller and
asks `allows`. Built-in names (`as_str`, `FromStr`) are what a command-line
flag accepts.

```rust
use openmassspec_core::{
    Check, Leniency, SpectrumRecord, Warning, WarningKind, WriteOptions,
};

// A built-in check, from the name a `--allow` flag carries.
let check: Check = "unknown_scan_mode".parse().unwrap();
let leniency = Leniency::strict().with(check);
let options = WriteOptions::new().with_leniency(leniency.clone());
assert!(options.leniency.allows(&Check::UnknownScanMode));

// A check that belongs to the reader.
let mine = Check::other("myreader.calibration_range");
assert!(!leniency.allows(&mine));
let leniency = leniency.with(mine.clone());
if leniency.allows(&mine) {
    let rec = SpectrumRecord::new(0, 1, "scan=1").with_warning(Warning::new(
        WarningKind::Other("myreader.calibration_range".into()),
        "m/z calibration is outside its stated range",
    ));
    assert_eq!(rec.warnings.len(), 1);
}
```

`Centroided::with_leniency` takes the same `Leniency`; with
`Check::UnknownScanMode` allowed it passes a spectrum of unknown scan mode
through unchanged and adds a `centroid_not_applied` warning.

#### A format that never states a value

Two checks guard values that some formats never state in any file: the
collision-energy unit (an Agilent Q-TOF scan record is a bare number) and
the scan mode (a SCIEX `.wiff` stream does not say profile or centroid).
A strict default would refuse every such file, so the reader declares the
fact once, in `RunMetadata::unstated`:

- `RunMetadata::with_unstated(check)` adds one check; `with_unstated_checks`
  sets the whole list; the field is `unstated: Vec<Check>`.
- Only `Check::UnknownScanMode` and `Check::UnknownCollisionEnergyUnit` are
  honoured. Any other entry is ignored.
- For a declared check the writer and `Centroided` behave as if the caller
  had allowed it: the value is written as in the lenient row of the table
  above and each affected record carries its warning
  (`unknown_scan_mode`, `unknown_collision_energy_unit`).
- The run also gets one `WarningKind::UnstatedBySource` warning (code
  `unstated_by_source`) per declared check that at least one record needed.
  It is written as a run-level `openmassspec: warning unstated_by_source`
  userParam and counted in `WriteSummary::warning_counts`, so
  `has_warnings` is true. A declared source whose records all state the
  value gets no warning.
- Declare it only when the format cannot state the value. A reader that
  reads the value from the file, or fails to for some files, declares
  nothing: a missing value then still fails loudly, which is how a decode
  bug shows up. The caller cannot switch a declaration off.

```rust
use openmassspec_core::{
    write_mzml, Check, CvTerm, PrecursorInfo, RunMetadata, ScanMode, SpectrumIter,
    SpectrumRecord, SpectrumSource,
};

struct QToF(Vec<SpectrumRecord>);

impl SpectrumSource for QToF {
    fn run_metadata(&self) -> RunMetadata {
        RunMetadata::new(
            "run.d",
            CvTerm::new("MS:1000563", "Thermo RAW format"),
            CvTerm::new("MS:1000768", "Thermo nativeID format"),
            CvTerm::new("MS:1000483", "Thermo Fisher Scientific instrument model"),
            "my-reader",
            "1.0.0",
        )
        // This format gives collision energies as bare numbers.
        .with_unstated(Check::UnknownCollisionEnergyUnit)
    }

    fn iter_spectra<'a>(&'a mut self) -> SpectrumIter<'a> {
        Box::new(self.0.clone().into_iter().map(Ok))
    }
}

let ms2 = SpectrumRecord::new(0, 1, "scan=1")
    .with_ms_level(2)
    .with_scan_mode(ScanMode::Centroid)
    .with_retention_time_sec(1.0)
    .with_precursor(
        PrecursorInfo::default()
            .with_selected_mz(500.0)
            .with_collision_energy(30.0, None),
    )
    .with_peaks(vec![100.0], vec![1.0]);

// Default options, no `Leniency`: the file converts, with warnings.
let summary = write_mzml(&mut QToF(vec![ms2]), &mut Vec::new()).unwrap();
assert_eq!(summary.warning_counts["unknown_collision_energy_unit"], 1);
assert_eq!(summary.warning_counts["unstated_by_source"], 1);
assert!(summary.has_warnings());
```

### 1.7 Removed `Copy` and `Eq` impls

`WriteSummary` is no longer `Copy` or `Eq`, and `WriteOptions` is no longer
`Eq`. Clone where you copied.

Before:

```rust,ignore
let summary = write_mzml(&mut src, &mut out)?;
log(summary);
if summary.skipped() > 0 { /* ... */ }   // used after the move
```

After:

```rust
use openmassspec_core::WriteSummary;

fn log(summary: WriteSummary) {
    let _ = summary;
}

let summary = WriteSummary::default();
log(summary.clone());
if summary.skipped() > 0 { /* ... */ }
```

### 1.8 `MsPower` is removed

Nothing in a record used it: `SpectrumRecord::ms_level` is a `u32`, and 4.0
rejects level 0. Use the number.

Before:

```rust,ignore
let level = MsPower::Ms2.ms_level();
```

After:

```rust
use openmassspec_core::SpectrumRecord;

let level: u32 = 2;
let rec = SpectrumRecord::new(0, 1, "scan=1").with_ms_level(level);
assert_eq!(rec.ms_level, 2);
```

### 1.9 Errors, and the Arrow module paths

`Error::Format` and `Error::Conformance` are removed. Wrap your own error
with `Error::vendor`, or report a problem with one record through Core with
`Error::record` and a namespaced code:

```rust
use openmassspec_core::{Error, RecordProblem};

#[derive(Debug)]
struct Truncated;
impl std::fmt::Display for Truncated {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("scan is truncated")
    }
}
impl std::error::Error for Truncated {}

let whole_file: Error = Error::vendor(Truncated);
// `Error::Vendor` exposes your error through `source()`, so callers can walk
// the chain and downcast.
assert!(std::error::Error::source(&whole_file).unwrap().is::<Truncated>());
let one_record: Error = Error::record(
    "scan=7",
    RecordProblem::other("myreader.bad_header", "header checksum mismatch"),
);
```

`SpectrumRecord::validate` and `ChromatogramRecord::validate` run the checks
the writer runs, so a reader's tests can assert that every record it yields
is writable. `conformance::assert_source_invariants` and
`assert_iter_invariants` are stricter than in 3.0.0; a reader's conformance
test can fail on data it passed before:

- it runs `SpectrumRecord::validate` on every record
  (`ConformanceError::InvalidRecord`);
- it rejects a repeated native id (`ConformanceError::DuplicateNativeId`);
- it requires a retention time on every spectrum or on none. For
  `assert_source_invariants` the run's `time_axis` decides (`Retention`: every
  spectrum; `Absent`: none); `assert_iter_invariants` has no run, so the first
  spectrum decides (`InvalidRecord` with `MissingRetentionTime` or
  `UnexpectedRetentionTime`);
- it rejects a declared base peak or m/z range on an empty spectrum, which
  it used to accept as 0.0 (`ConformanceError::DeclaredOnEmptySpectrum`), and
  an empty spectrum without a `WarningKind::EmptySpectrum` warning
  (`ConformanceError::EmptySpectrum`); see 1.11.

If you use `openmassspec-arrow`, everything moves into a module named for
the arrow-rs release it is built against (feature `arrow-60`, on by
default). Only version-independent items stay at the crate root:
`SCHEMA_VERSION`, `SCHEMA_VERSION_KEY`, `COLUMNS`,
`check_schema_version_value` and `SchemaVersionError`. A newer arrow-rs
major arrives as a new module beside `arrow_60`, so these paths stay valid.
`push` validates the record and returns a `Result`, and `BatchError` (one
per module) replaces `ArrowError` as the error type. `push` rejects what
`SpectrumRecord::validate` rejects and a native id the builder has already
seen, and returns `BatchError::ColumnRange` for an `index` that does not fit
`u32` or an `ms_level` that does not fit `u8`, which 3.0.0 truncated with
`as`.

Before:

```rust,ignore
use openmassspec_arrow::{spectrum_record_schema, SpectrumBatchBuilder};

let mut builder = SpectrumBatchBuilder::new(None);
for rec in &records {
    builder.push(rec);
}
let batch = builder.finish()?;   // Result<RecordBatch, ArrowError>
```

After:

```rust
use openmassspec_arrow::arrow_60::{spectrum_record_schema, BatchError, SpectrumBatchBuilder};
use openmassspec_core::SpectrumRecord;

fn rows(records: &[SpectrumRecord]) -> Result<usize, BatchError> {
    let mut builder = SpectrumBatchBuilder::new(None);
    for rec in records {
        builder.push(rec)?; // validates the record and the native id
    }
    let batch = builder.finish()?;
    assert_eq!(batch.schema(), spectrum_record_schema());
    Ok(batch.num_rows())
}

let record = SpectrumRecord::new(0, 1, "scan=1").with_peaks(vec![100.0], vec![5.0]);
assert_eq!(rows(&[record]).unwrap(), 1);
```

### 1.10 Output that changes with no compile error

Some changes show up only in the mzML. The complete list is in 2.3. Two of
them belong to the reader:

- `Activation::CID` is always written as collision-induced dissociation
  (`MS:1000133`). 3.0.0 wrote beam-type CID (`MS:1000422`) when the analyzer
  was `FTMS`, so a reader that sent `CID` for Orbitrap beam-type scans now
  gets resonant CID in the file. A reader that knows a scan is beam-type
  sends `Activation::HCD`; one that knows it is resonant keeps `CID`.
  `PrecursorInfo::analyzer` no longer chooses the term.
- `<scanWindowList>` comes only from `scan_windows` (1.5), and a precursor
  whose activation is `Activation::PD` or unset gets an
  `openmassspec: activation` userParam (`pd` or `unknown`).

```rust
use openmassspec_core::{Activation, Analyzer, PrecursorInfo};

// An Orbitrap beam-type scan: HCD, not CID.
let beam = PrecursorInfo::default()
    .with_activation(Activation::HCD)
    .with_analyzer(Analyzer::FTMS);
assert_eq!(beam.activation, Some(Activation::HCD));
```

### 1.11 Empty spectra

A spectrum with no peaks needs two things from the reader, and the
conformance harness checks both:

- It carries `Warning::new(WarningKind::EmptySpectrum, message)` (code
  `empty_spectrum`), so a consumer finds every empty spectrum by one code
  whichever reader produced it. A reader-specific code such as
  `myreader.empty_spectrum` does not count; add it as well if you want it.
  The flag reaches the mzML as a warning userParam and the Arrow batch in
  the `warnings` column, and it makes `WriteSummary::has_warnings` true.
  Do not treat that as a failure: an empty scan is a fact about the run.
- It declares no base peak and no m/z range: `base_peak_mz`,
  `base_peak_intensity`, `low_mz` and `high_mz` stay `None`. A declared
  `0.0` is rejected by the harness
  (`ConformanceError::DeclaredOnEmptySpectrum`). `effective_base_peak` and
  `effective_mz_range` return `None` for an empty spectrum even when those
  fields are set, so the writer never prints them. `total_ion_current` may
  be declared as `0.0`.

The writer does not require the flag; the harness does, so a reader that
forgets it fails in its own CI.

```rust
use openmassspec_core::conformance::assert_iter_invariants;
use openmassspec_core::{SpectrumRecord, Warning, WarningKind};

let empty = SpectrumRecord::new(0, 1, "scan=1")
    .with_retention_time_sec(1.0)
    .with_warning(Warning::new(WarningKind::EmptySpectrum, "the scan has no peaks"));
assert_eq!(empty.effective_base_peak(), None);
assert_eq!(empty.effective_mz_range(), None);
assert_eq!(assert_iter_invariants(vec![empty.clone()]).unwrap(), 1);

// Without the flag the harness refuses it.
let mut unflagged = empty.clone();
unflagged.warnings.clear();
assert!(assert_iter_invariants(vec![unflagged]).is_err());

// So does a declared base peak, even 0.0.
let declared = empty.with_base_peak_intensity(0.0);
assert!(assert_iter_invariants(vec![declared]).is_err());
```

## 2. Users of records

You use the records a reader returns, or the writer, from Rust.

### 2.1 Field changes

| 3.0.0                                        | 4.0.0                                                              |
| -------------------------------------------- | ------------------------------------------------------------------ |
| `PrecursorInfo::collision_energy: Option<f64>` | `PrecursorInfo::collision_energies: Vec<CollisionEnergy>`         |
| `PrecursorInfo::ce_is_nce: bool`             | `CollisionEnergy::unit: Option<CollisionEnergyUnit>` on each entry |
| `SpectrumRecord::retention_time_sec: f64`    | `Option<f64>`                                                      |
| (none)                                       | `SpectrumRecord::scan_windows: Vec<ScanWindow>`                    |
| (none)                                       | `warnings: Vec<Warning>` on `SpectrumRecord`, `ChromatogramRecord` and `RunMetadata` |
| (none)                                       | `RunMetadata::time_axis: TimeAxis`                                 |
| `MsPower`                                    | removed; `ms_level` is a `u32`                                     |

Before:

```rust,ignore
let nce = if p.ce_is_nce { p.collision_energy } else { None };
let minutes = rec.retention_time_sec / 60.0;
```

After:

```rust
use openmassspec_core::{CollisionEnergyUnit, PrecursorInfo, SpectrumRecord};

fn nce(p: &PrecursorInfo) -> Option<f64> {
    p.collision_energies
        .iter()
        .find(|ce| ce.unit == Some(CollisionEnergyUnit::Nce))
        .map(|ce| ce.value)
}

fn minutes(rec: &SpectrumRecord) -> Option<f64> {
    rec.retention_time_sec.map(|seconds| seconds / 60.0)
}

let p = PrecursorInfo::default().with_collision_energy(30.0, CollisionEnergyUnit::Nce);
assert_eq!(nce(&p), Some(30.0));
assert_eq!(minutes(&SpectrumRecord::new(0, 1, "scan=1")), None);
```

An energy with `unit: None` is a value the reader could not place. Treat it
as unknown instead of as eV.

### 2.2 Serialized labels

Serde writes and reads the `as_str` labels of `Polarity`, `ScanMode`,
`Analyzer`, `Activation`, `MobilityArrayKind` and `CollisionEnergyUnit`, not
the Rust variant names. `Analyzer::TOFMS` is `"tof"` and `Polarity::Positive`
is `"positive"`. Each type has `FromStr` for the same labels. JSON written by
Core 3 (`"Positive"`, `"TOFMS"`) does not deserialize, and neither does a
precursor with a `collision_energy` or `ce_is_nce` key: the error names the
replacement. Convert stored records, or read them again from the source file.

```rust
use openmassspec_core::{Analyzer, Polarity};

assert_eq!(Analyzer::TOFMS.as_str(), "tof");
assert_eq!("positive".parse::<Polarity>().unwrap(), Polarity::Positive);
assert!("TOFMS".parse::<Analyzer>().is_err());
```

### 2.3 Conversions are stricter

Records that Core 3 wrote with a silent default are errors now (unless the
source declares that its format never states the value, 1.6): an
unknown scan mode, a collision energy with no unit, a missing retention
time, ms level 0, unsorted or non-finite m/z, mismatched arrays, a repeated
native id, or a chromatogram with mismatched arrays, non-finite values or a
repeated id (none of the array and numeric problems has a `Check`, so a
caller cannot allow them), and a precursor reference to a spectrum that has not been
written yet (the writer resolves `precursor_native_id` against spectra it
has already written, so stream order matters; see 1.6). The record checks
and what each asks of the reader are in 1.5. The "Behavior changes" list in
[CHANGELOG.md](CHANGELOG.md) has all of them.

A source whose `spectrum_count_hint` is larger than the number of spectra it
yields, with no `Err` for the difference, is `Error::SpectrumCountShortfall`
unless `Check::SpectrumCountShortfall` is allowed (1.6), so a reader's hint
must not exceed what it yields, or be `None`.

`WriteSummary` says what happened. With `RecordErrorPolicy::Skip` a
rejected record is left out and counted, and its cause is kept in `issues`:

```rust
use openmassspec_core::{
    write_mzml_with_options, RecordErrorPolicy, SpectrumSource, WriteOptions,
};

fn convert<S: SpectrumSource>(source: &mut S, out: &mut Vec<u8>) -> openmassspec_core::Result<()> {
    let options = WriteOptions::new().with_record_error_policy(RecordErrorPolicy::Skip);
    let summary = write_mzml_with_options(source, out, &options)?;
    for issue in &summary.issues {
        eprintln!("{:?} {}: {}", issue.kind, issue.position, issue.cause);
    }
    if summary.skipped() > 0 {
        eprintln!("{} records skipped", summary.skipped());
    }
    Ok(())
}
```

To accept a specific kind of record instead of skipping, allow its `Check`
(1.6). `WriteSummary::has_warnings` is true when anything was skipped,
flagged or left unreconciled; `warning_counts` counts the warnings per code.

Output that changes without an error, for anyone who compares mzML from 3.x
and 4.x:

- `Activation::CID` is always `MS:1000133`. 3.0.0 wrote `MS:1000422`
  (beam-type) for it when the analyzer was `FTMS`; readers that know a scan
  is beam-type now send `HCD`, so the same file can differ in the activation
  term.
- `<scanWindowList>` is written only from `scan_windows`. It used to carry
  the observed m/z range. Empty spectra omit the base-peak and observed-range
  terms instead of writing 0, even when the record declares them, and a
  spectrum with no retention time has no scan start time.
- A spectrum whose analyzer is not in `RunMetadata::analyzers` gets its own
  `instrumentConfiguration` (it used to inherit `IC1`). `<fileContent>`
  declares only the spectrum kinds written.
- Precursors with `Activation::PD` or no activation carry an
  `openmassspec: activation` userParam, and stepped collision energies are
  all written.
- Warnings are written as `openmassspec: warning CODE` userParams.
- The writer writes the document head after it has read and checked all
  spectra. A spectrum that fails (a source `Err`, a rejected record, the
  spectrum-count check, the spool file cannot be created or written)
  therefore leaves `out` empty, for a `Vec` and for a `File`; 3.0.0 left a
  truncated document. A failure after the head is written still leaves an
  incomplete document, with no `</run>`: a chromatogram error leaves the head
  and the complete `<spectrumList>`, an I/O error from `out` leaves what it
  accepted, and a failure reading the spool back (flush, read, byte count,
  removal) leaves the head and an open `<spectrumList>`, with none, some or all
  spectra. Write to a temporary file and rename it on success if you need all
  or nothing.
- `SpectrumRecord::effective_base_peak` reports a NaN intensity instead of
  skipping it when it comes first.

### 2.4 Arrow batches

The module paths and error type are in 1.9. The columns change as follows:

| Column                                | 3.0.0              | 4.0.0                                                    |
| ------------------------------------- | ------------------ | -------------------------------------------------------- |
| `retention_time_sec`                  | `Float64`, not null | `Float64`, nullable                                      |
| `precursor_ce_is_nce` (24th column, index 23) | `UInt8`    | replaced in place by `precursor_collision_energy_unit`   |
| `precursor_collision_energy_unit` (index 23)  | -          | `Utf8` (`"ev"`, `"nce"`), null when unknown; follows `precursor_collision_energy` |
| `precursor_collision_energies`        | -                  | `LargeList<Struct{value: Float64, unit: Utf8}>`, every energy |
| `scan_windows`                        | -                  | `LargeList<Struct{lower: Float64, upper: Float64}>`      |
| `warnings`                            | -                  | `LargeList<Struct{code: Utf8, message: Utf8}>` (spectrum warnings) |

A consumer that finds columns by position instead of by name breaks at index
23, where a `UInt8` flag became a `Utf8` label. The three new list columns
are appended at the end. `precursor_collision_energy` and its unit column
hold the first energy;
`precursor_collision_energies` holds all of them. The schema metadata
carries `openmassspec.schema_version` (currently 1). A reader accepts that
version or older and rejects a newer one. Arrow data written by the 3.x
crate has no version key and is rejected with an error that says so: write
it again from the source file.

## 3. Python users

OpenTFRaw, OpenWRaw and OpenTimsTDF build their record dicts from Core's
serialized form, so the changes in 2.1 and 2.2 show up as dict keys, and
the table below describes those three. OpenARaw, OpenSXRaw, OpenSZRaw and
the `openmassspec` package build their dicts by hand and decide their own
keys and labels in their own releases; their CHANGELOGs and migration notes
document them.

| Before                                                  | After                                                                  |
| ------------------------------------------------------- | ---------------------------------------------------------------------- |
| `rec["precursor"]["collision_energy"]` (float or None)  | `rec["precursor"]["collision_energies"]`, a list of `{"value", "unit"}` |
| `rec["precursor"]["ce_is_nce"]` (bool)                  | the `"unit"` of each entry: `"ev"`, `"nce"` or `None`                  |
| `rec["retention_time_sec"]` is always a float           | may be `None`                                                          |
| `analyzer` is `"tofms"` (also `precursor["analyzer"]` and the run's `analyzers`) | `"tof"`; the other labels are unchanged |
| no `scan_windows`                                       | `rec["scan_windows"]`, a list of `{"lower", "upper"}`                  |
| no `warnings`                                           | `rec["warnings"]`, a list of `{"kind", "message"}`                     |
| run metadata without `time_axis`                        | `"retention_time"` or `"absent"`                                       |
| no `unstated` in run metadata                           | `unstated`, a list of check names (`"unknown_scan_mode"`, `"unknown_collision_energy_unit"`); left out when empty |

`collision_energies`, `scan_windows` and `warnings` are left out of the dict
when empty, so read them with `.get(key, [])`. An empty spectrum has a
warning with kind `"empty_spectrum"` and no base peak or m/z range. An MS1 spectrum has
`precursor` set to `None`.

Before:

```python
nce = None
precursor = rec["precursor"]
if precursor is not None and precursor["ce_is_nce"]:
    nce = precursor["collision_energy"]

minutes = rec["retention_time_sec"] / 60.0
is_tof = rec["analyzer"] == "tofms"
```

After:

```python
precursor = rec["precursor"]
energies = precursor.get("collision_energies", []) if precursor is not None else []
nce = [e["value"] for e in energies if e["unit"] == "nce"]
unknown_unit = [e["value"] for e in energies if e["unit"] is None]

rt = rec["retention_time_sec"]
minutes = rt / 60.0 if rt is not None else None
is_tof = rec["analyzer"] == "tof"
```

Conversions to mzML and the command line fail on the records listed in 2.3,
where Core 3 wrote them. A reader that exposes the checks lets the caller
opt in by name (`unknown_scan_mode`, `dangling_precursor_ref`,
`unknown_collision_energy_unit`, `spectrum_count_shortfall`); the reader's
documentation says how.

The six reader packages do not return Arrow. Arrow output of the bridge
follows 2.4.
