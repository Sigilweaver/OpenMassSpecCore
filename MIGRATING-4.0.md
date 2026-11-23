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
you. Every `rust` block that follows a "After" heading is compiled as a doc
test in CI; blocks marked "Before" show 3.0.0 code and are not compiled.

Three changes decide most of the work:

- Collision energy is a list, and each entry names its unit or says it is
  unknown. Core no longer assumes eV.
- Retention time is optional, and a spectrum with no scan mode is an error.
  A record that Core 3 wrote with a silent default is now rejected.
- Iterators are `Send`, so a reader generic over `Read + Seek` adds `Send`
  to its bounds. That is a public signature change in the reader, so the
  reader needs a major release too.

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
The iterator borrows `&mut self`, so everything the reader holds must be
`Send`. A reader generic over its byte source adds `Send` to the bound, on
the struct and on the impl. A field that is not `Send` (`Rc`, `RefCell`,
a raw pointer) is replaced by `Arc`, `Mutex` or a plain owned value.

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
as a unit-less `userParam` and adds a warning to the record.

### 1.4 Retention time is optional

`SpectrumRecord::retention_time_sec` is `Option<f64>`, and
`SpectrumRecord::new` leaves it `None` (it used to set `0.0`). The setter
still takes seconds as an `f64` and also accepts an `Option<f64>`.

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
with `Check::UnknownScanMode` (1.6). `Centroided` picks peaks only from
spectra declared `ScanMode::Profile`.

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

### 1.6 Leniency, `Check` and warnings

By default the writer and `Centroided` reject what they cannot represent
faithfully. A `Leniency` names the `Check`s that may continue instead, and
each affected record gets a `Warning`:

| Check                          | Strict mode                                | Lenient mode                                    |
| ------------------------------ | ------------------------------------------ | ----------------------------------------------- |
| `UnknownScanMode`              | record error                               | written with no profile or centroid term        |
| `DanglingPrecursorRef`         | record error                               | the `spectrumRef` is dropped                    |
| `UnknownCollisionEnergyUnit`   | record error                               | value written as a unit-less `userParam`        |
| `SpectrumCountShortfall`       | `Error::SpectrumCountShortfall`            | run-level warning and `WriteSummary::spectrum_count_mismatch` |

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

`Centroided::with_leniency` takes the same `Leniency`.

### 1.7 Removed `Copy` and `Eq` impls

`WriteSummary` is no longer `Copy` or `Eq`, and `WriteOptions` is no longer
`Eq`. The types added in 4.0 (`Check`, `Leniency`, `ScanWindow`,
`CollisionEnergy`, `CountMismatch`) are not `Copy` either. Clone where you
copied, and pass `&Check` to `Leniency::allows`.

Before:

```rust,ignore
let summary = write_mzml(&mut src, &mut out)?;
log(summary);
if summary.skipped() > 0 { /* ... */ }   // used after the move
```

After:

```rust
use openmassspec_core::{Check, Leniency, WriteSummary};

fn log(summary: WriteSummary) {
    let _ = summary;
}

let summary = WriteSummary::default();
log(summary.clone());
if summary.skipped() > 0 { /* ... */ }

let leniency = Leniency::strict().with(Check::DanglingPrecursorRef);
assert!(leniency.allows(&Check::DanglingPrecursorRef));
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
let one_record: Error = Error::record(
    "scan=7",
    RecordProblem::other("myreader.bad_header", "header checksum mismatch"),
);
```

`SpectrumRecord::validate` and `ChromatogramRecord::validate` run the checks
the writer runs, so a reader's tests can assert that every record it yields
is writable. `conformance::assert_source_invariants` applies them too.

If you use `openmassspec-arrow`, everything moves into a module named for
the arrow-rs release it is built against (feature `arrow-60`, on by
default). Only version-independent items stay at the crate root:
`SCHEMA_VERSION`, `SCHEMA_VERSION_KEY`, `COLUMNS`,
`check_schema_version_value` and `SchemaVersionError`. A newer arrow-rs
major arrives as a new module beside `arrow_60`, so these paths stay valid.
`push` validates the record and returns a `Result`, and `BatchError` (one
per module) replaces `ArrowError` as the error type.

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

Records that Core 3 wrote with a silent default are errors now: an
unknown scan mode, a collision energy with no unit, a missing retention
time, ms level 0, unsorted or non-finite m/z, mismatched arrays, a repeated
native id, and a precursor reference to a spectrum that is not in the
output. The "Behavior changes" list in [CHANGELOG.md](CHANGELOG.md) has all
of them.

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
    if summary.has_warnings() {
        eprintln!("{} records skipped", summary.skipped());
    }
    Ok(())
}
```

To accept a specific kind of record instead of skipping, allow its `Check`
(1.6).

### 2.4 Arrow batches

The module paths and error type are in 1.9. The columns change as follows:

| Column                                | 3.0.0              | 4.0.0                                                    |
| ------------------------------------- | ------------------ | -------------------------------------------------------- |
| `retention_time_sec`                  | `Float64`, not null | `Float64`, nullable                                      |
| `precursor_ce_is_nce`                 | `UInt8`            | removed                                                  |
| `precursor_collision_energy_unit`     | -                  | `Utf8` (`"ev"`, `"nce"`), null when unknown; follows `precursor_collision_energy` |
| `precursor_collision_energies`        | -                  | `LargeList<Struct{value: Float64, unit: Utf8}>`, every energy |
| `scan_windows`                        | -                  | `LargeList<Struct{lower: Float64, upper: Float64}>`      |
| `warnings`                            | -                  | `LargeList<Struct{code: Utf8, message: Utf8}>` (spectrum warnings) |

`precursor_collision_energy` and its unit column hold the first energy;
`precursor_collision_energies` holds all of them. The schema metadata
carries `openmassspec.schema_version` (currently 1). A reader accepts that
version or older and rejects a newer one. Arrow data written by the 3.x
crate has no version key and is rejected with an error that says so: write
it again from the source file.

## 3. Python users

OpenTFRaw, OpenWRaw and OpenTimsTDF build their record dicts from Core's
serialized form, so the changes in 2.1 and 2.2 show up as dict keys.
OpenARaw, OpenSXRaw and OpenSZRaw build theirs by hand. The table is the
shape Core produces; the CHANGELOG of each reader names the release that
returns it.

| Before                                                  | After                                                                  |
| ------------------------------------------------------- | ---------------------------------------------------------------------- |
| `rec["precursor"]["collision_energy"]` (float or None)  | `rec["precursor"]["collision_energies"]`, a list of `{"value", "unit"}` |
| `rec["precursor"]["ce_is_nce"]` (bool)                  | the `"unit"` of each entry: `"ev"`, `"nce"` or `None`                  |
| `rec["retention_time_sec"]` is always a float           | may be `None`                                                          |
| `analyzer` is `"tofms"` (also `precursor["analyzer"]` and the run's `analyzers`) | `"tof"`; the other labels are unchanged |
| no `scan_windows`                                       | `rec["scan_windows"]`, a list of `{"lower", "upper"}`                  |
| no `warnings`                                           | `rec["warnings"]`, a list of `{"kind", "message"}`                     |
| run metadata without `time_axis`                        | `"retention_time"` or `"absent"`                                       |

`collision_energies`, `scan_windows` and `warnings` are left out of the dict
when empty, so read them with `.get(key, [])`. An MS1 spectrum has
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

Arrow output follows 2.4. The six reader packages do not return Arrow.
`read_arrow` in the `openmassspec` package returns what the bridge builds,
so a build against Core 4 has the new columns and the schema version, and
rejects Arrow files written with the 3.x bridge.
