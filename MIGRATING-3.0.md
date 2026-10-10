# Migrating a reader to openmassspec-core 3.0

Mechanical steps for a vendor reader crate that implements `SpectrumSource`
and builds `SpectrumRecord` / `RunMetadata` / `ChromatogramRecord` values.
Work through them in order; `cargo check` after each step points at the
remaining sites.

## 1. Bump the dependency

```toml
[dependencies]
openmassspec-core = "3.0.0"
```

If you used the `arrow` feature, drop it and add the new crate instead:

```toml
openmassspec-arrow = "3.0.0"
```

and replace `openmassspec_core::arrow::` with `openmassspec_arrow::`. The
function and type names are unchanged.

## 2. Yield `Result` items from `iter_spectra` / `iter_chromatograms`

The trait now returns `SpectrumIter<'a>` / `ChromatogramIter<'a>`, aliases
for `Box<dyn Iterator<Item = openmassspec_core::Result<_>> + 'a>`.

Wrap successful records in `Ok` and turn each decode failure into an `Err`
instead of skipping it. `Error::vendor` boxes any
`std::error::Error + Send + Sync + 'static`.

Before:

```rust
fn iter_spectra<'s>(&'s mut self) -> Box<dyn Iterator<Item = msc::SpectrumRecord> + 's> {
    Box::new(std::iter::from_fn(move || {
        while idx < n {
            let cur = idx;
            idx += 1;
            if let Some(parts) = extract_parts(raw, source, cur) {
                return Some(to_msc_record(parts));
            }
            // unreadable scan: silently skipped
        }
        None
    }))
}
```

After:

```rust
fn iter_spectra<'s>(&'s mut self) -> msc::SpectrumIter<'s> {
    Box::new(std::iter::from_fn(move || {
        if idx >= n {
            return None;
        }
        let cur = idx;
        idx += 1;
        Some(
            try_extract_parts(raw, source, cur)
                .map(to_msc_record)
                .map_err(msc::Error::vendor),
        )
    }))
}
```

If a scan is legitimately absent (not an error, e.g. a controller with no
data for that index), keep skipping it. Only failures become `Err`.

If your internal decoder returns `Option` and throws the error away, change
it to return `Result<_, YourError>` first. A quick bridge while you do that:

```rust
.ok_or_else(|| msc::Error::vendor(YourError::ScanUnreadable(cur)))
```

For chromatograms built eagerly into a `Vec`:

```rust
fn iter_chromatograms<'s>(&'s mut self) -> msc::ChromatogramIter<'s> {
    Box::new(build_chromatograms(self.raw).into_iter().map(Ok))
}
```

Optional: implement `From<YourError> for openmassspec_core::Error` (allowed
by the orphan rules because `YourError` is local) so `?` works inside
helpers that return `openmassspec_core::Result`:

```rust
impl From<YourError> for openmassspec_core::Error {
    fn from(e: YourError) -> Self {
        openmassspec_core::Error::vendor(e)
    }
}
```

## 3. Build records with constructors and setters

`SpectrumRecord`, `PrecursorInfo`, `ChromatogramRecord`, `RunMetadata`, and
`CvTerm` are `#[non_exhaustive]`: struct literals no longer compile outside
the core crate. Fields stay public, so reading and assigning them works.

Every field has a `with_<field>` setter. Setters for `Option<T>` fields
accept either `T` or `Option<T>`, so most struct-literal lines become one
setter call with the same right-hand side.

`SpectrumRecord` before:

```rust
msc::SpectrumRecord {
    index: meta.index,
    scan_number: meta.scan_number,
    native_id: native_id_for(meta.scan_number),
    ms_level: meta.ms_level,
    polarity: convert_polarity(meta.polarity),
    retention_time_sec: meta.retention_time_min * 60.0,
    total_ion_current: Some(meta.total_ion_current),
    faims_cv: meta.faims_cv,
    inv_mobility: None,
    precursor,
    mz,
    intensity,
    inv_mobility_per_peak: None,
    extra,
    // ...every other field...
}
```

After (fields that were `None` / empty can be dropped; `new` sets them):

```rust
msc::SpectrumRecord::new(meta.index, meta.scan_number, native_id_for(meta.scan_number))
    .with_ms_level(meta.ms_level)
    .with_polarity(convert_polarity(meta.polarity))
    .with_retention_time_sec(meta.retention_time_min * 60.0)
    .with_total_ion_current(meta.total_ion_current)
    .with_faims_cv(meta.faims_cv)
    .with_precursor(precursor)
    .with_peaks(mz, intensity)
    .with_extra_map(extra)
```

`SpectrumRecord::new` sets `ms_level` to 1, retention time to 0, empty
peaks, and every optional field to `None`.

`PrecursorInfo` before / after:

```rust
msc::PrecursorInfo { target_mz: p.target_mz, charge: p.charge, ce_is_nce: p.ce_is_nce, intensity: None, ccs: None, /* ... */ }

msc::PrecursorInfo::default()
    .with_target_mz(p.target_mz)
    .with_charge(p.charge)
    .with_ce_is_nce(p.ce_is_nce)
```

`ChromatogramRecord` before / after:

```rust
msc::ChromatogramRecord { index: 0, id: "TIC".into(), chromatogram_type: Some(term), precursor_mz: None, product_mz: None, time_sec, intensity }

msc::ChromatogramRecord::new(0, "TIC", time_sec, intensity).with_chromatogram_type(term)
```

`RunMetadata` before / after:

```rust
msc::RunMetadata { source_file_name, source_file_format, native_id_format, instrument, software_name, software_version, start_timestamp, analyzers, extra, /* None fields */ }

msc::RunMetadata::new(source_file_name, source_file_format, native_id_format, instrument, software_name, software_version)
    .with_start_timestamp(start_timestamp)
    .with_analyzers(analyzers)
    .with_extra_map(extra)
```

Test code that used `SpectrumRecord { field: x, ..base }` becomes
`let mut r = base; r.field = x;` or `base.with_field(x)`.

## 4. Add wildcard arms when matching core enums

`Polarity`, `ScanMode`, `Analyzer`, `MsPower`, `Activation`, `Compression`,
`MobilityArrayKind`, `Error`, and `conformance::ConformanceError` are
`#[non_exhaustive]`. Conversions *into* core enums (`match vendor_enum {
.. => msc::Analyzer::FTMS }`) need no change. A `match` *on* a core enum
needs a `_ =>` arm; decide explicitly what an unknown variant means for your
code (usually an error or `None`, not a silent default).

## 5. `CvTerm` fields are `Cow<'static, str>`

`CvTerm::new("MS:1000563", "Thermo RAW format")` still works with string
literals and also accepts `String`. Code that read `term.accession` as a
`&'static str` should use `&term.accession` / `term.accession.as_ref()`.
Comparisons such as `term.accession == "MS:1000563"` still compile.

For terms the core writer already knows, prefer the typed constants:
`openmassspec_core::cv::TOTAL_ION_CURRENT_CHROMATOGRAM.into()`.

## 6. Writer calls return `Result<WriteSummary, openmassspec_core::Error>`

`write_mzml` / `write_indexed_mzml` now return
`openmassspec_core::Result<WriteSummary>` (was `std::io::Result<()>`).
`*_with_compression` is replaced by `*_with_options`:

```rust
// before
openmassspec_core::write_mzml_with_compression(&mut src, out, Compression::NoCompression)?;

// after
let opts = WriteOptions::new().with_compression(Compression::NoCompression);
openmassspec_core::write_mzml_with_options(&mut src, out, &opts)?;
```

If your crate's error type wrapped `std::io::Error` from the writer, add a
variant (or `From` impl) for `openmassspec_core::Error`:

```rust
#[error(transparent)]
Core(#[from] openmassspec_core::Error),
```

The writer stops at the first `Err` record. A CLI that wants best-effort
output can opt in and must surface the count:

```rust
let opts = WriteOptions::new().with_record_error_policy(RecordErrorPolicy::Skip);
let summary = openmassspec_core::write_mzml_with_options(&mut src, out, &opts)?;
if summary.skipped() > 0 {
    eprintln!("warning: {} records could not be decoded and were skipped", summary.skipped());
}
```

## 7. `additional_processing_steps` returns `Vec<CvTerm>`

Only adapters override this. `vec![("MS:1000035", "peak picking")]` becomes
`vec![openmassspec_core::cv::PEAK_PICKING.into()]`.

## 8. Consumers of `iter_spectra`

Code that iterates a source directly now gets `Result` items:

```rust
for rec in source.iter_spectra() {
    let rec = rec?;
    // ...
}
// or
let recs: Vec<_> = source.iter_spectra().collect::<openmassspec_core::Result<_>>()?;
```

`conformance::assert_source_invariants` reports an `Err` item as
`ConformanceError::Source { position, error }`.

## 9. Update docs

Remove any doc line that says unreadable scans are "skipped silently"; state
that they are reported as errors.
