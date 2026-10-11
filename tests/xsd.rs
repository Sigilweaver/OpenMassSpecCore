//! Validates writer output against the vendored PSI-MS mzML XSDs with
//! `xmllint`, covering shapes the CI `validate-mzml` example does not: a
//! digit-leading source file name (the `<run id>` is an `xs:ID`), a software
//! name that is not an NCName (the `<software id>` is an `xs:ID`), an MS2
//! precursor with unknown activation, an SRM chromatogram precursor (both of
//! which write an empty `<activation>`), and a document written with skipped
//! records.
//!
//! Skips with a message when `xmllint` is not installed.

use std::path::{Path, PathBuf};
use std::process::Command;

use openmassspec_core::{
    write_indexed_mzml, write_indexed_mzml_with_options, write_mzml, write_mzml_with_options,
    ChromatogramIter, ChromatogramRecord, CvTerm, Error, Polarity, PrecursorInfo,
    RecordErrorPolicy, RunMetadata, ScanMode, SpectrumIter, SpectrumRecord, SpectrumSource,
    WriteOptions,
};

struct Source {
    meta: RunMetadata,
    spectra: Vec<SpectrumRecord>,
    chroms: Vec<ChromatogramRecord>,
    /// Yield an `Err` between the spectra.
    with_error: bool,
    /// Yield an `Err` after the last spectrum and in place of the
    /// chromatogram.
    with_more_errors: bool,
    /// Return `None` from `spectrum_count_hint`.
    no_hint: bool,
}

fn spectrum(index: usize, ms_level: u32, precursor: Option<PrecursorInfo>) -> SpectrumRecord {
    SpectrumRecord::new(
        index,
        index as u32 + 1,
        format!("controllerType=0 controllerNumber=1 scan={}", index + 1),
    )
    .with_ms_level(ms_level)
    .with_polarity(Polarity::Positive)
    .with_scan_mode(ScanMode::Centroid)
    .with_retention_time_sec(index as f64)
    .with_precursor(precursor)
    .with_peaks(vec![100.0, 200.0], vec![1.0, 2.0])
}

impl Source {
    fn new(source_file_name: &str, software_name: &str) -> Self {
        Self {
            meta: RunMetadata::new(
                source_file_name,
                CvTerm::new("MS:1000563", "Thermo RAW format"),
                CvTerm::new("MS:1000768", "Thermo nativeID format"),
                CvTerm::new("MS:1001911", "Q Exactive"),
                software_name,
                "0.0.0",
            )
            .with_acquisition_software_name("Xcalibur".to_string())
            .with_acquisition_software_version("4.4".to_string()),
            spectra: vec![
                spectrum(0, 1, None),
                spectrum(
                    1,
                    2,
                    Some(
                        PrecursorInfo::default()
                            .with_selected_mz(150.0)
                            .with_collision_energy(30.0)
                            .with_ce_is_nce(true),
                    ),
                ),
            ],
            chroms: vec![ChromatogramRecord::new(
                0,
                "SRM SIC Q1=524.3 Q3=136.1",
                vec![0.0, 30.0],
                vec![5.0, 8.0],
            )
            .with_chromatogram_type(CvTerm::new(
                "MS:1001473",
                "selected reaction monitoring chromatogram",
            ))
            .with_precursor_mz(524.3)
            .with_product_mz(136.1)],
            with_error: false,
            with_more_errors: false,
            no_hint: false,
        }
    }
}

impl SpectrumSource for Source {
    fn run_metadata(&self) -> RunMetadata {
        self.meta.clone()
    }
    fn iter_spectra<'a>(&'a mut self) -> SpectrumIter<'a> {
        let error = self
            .with_error
            .then(|| Err(Error::Format("corrupt scan".into())));
        let trailing = self
            .with_more_errors
            .then(|| Err(Error::Format("corrupt last scan".into())));
        let mut spectra = self.spectra.iter().cloned().map(Ok);
        let first = spectra.next();
        Box::new(
            first
                .into_iter()
                .chain(error)
                .chain(spectra)
                .chain(trailing),
        )
    }
    fn iter_chromatograms<'a>(&'a mut self) -> ChromatogramIter<'a> {
        let error = self
            .with_more_errors
            .then(|| Err(Error::Format("corrupt chromatogram".into())));
        Box::new(
            error
                .into_iter()
                .chain(self.chroms.clone().into_iter().map(Ok)),
        )
    }
    fn spectrum_count_hint(&self) -> Option<usize> {
        // Counts the failing records too, as a reader that knows its scan
        // count but not which scans will fail would.
        let failing = usize::from(self.with_error) + usize::from(self.with_more_errors);
        (!self.no_hint).then_some(self.spectra.len() + failing)
    }
}

fn xmllint_available() -> bool {
    Command::new("xmllint")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
}

fn validate(doc: &[u8], schema: &str, name: &str) {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
    let path = dir.join(name);
    std::fs::write(&path, doc).unwrap();
    let schema = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("schema")
        .join(schema);
    let out = Command::new("xmllint")
        .arg("--noout")
        .arg("--schema")
        .arg(&schema)
        .arg(&path)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{name} failed XSD validation:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn awkward_file_and_software_names_validate_against_xsd() {
    if !xmllint_available() {
        eprintln!("skipping: xmllint not installed");
        return;
    }
    for (file, software) in [
        ("20260101_sample.raw", "xsd-test"),
        ("1 run (copy).raw", "1 reader & co"),
        ("same", "same"),
        ("acquisition_software", "acquisition_software"),
    ] {
        let mut src = Source::new(file, software);

        let mut plain = Vec::new();
        write_mzml(&mut src, &mut plain).unwrap();
        validate(&plain, "mzML1.1.0.xsd", "xsd_plain.mzML");

        let mut indexed = Vec::new();
        write_indexed_mzml(&mut src, &mut indexed).unwrap();
        validate(&indexed, "mzML1.1.2_idx.xsd", "xsd_indexed.mzML");
    }
}

/// The `count` attribute of `<{list}`, the `index` attributes of every
/// `<{elem}`, in document order.
fn count_and_indexes(doc: &str, list: &str, elem: &str) -> (Option<usize>, Vec<usize>) {
    fn attr(tag: &str, name: &str) -> usize {
        let key = format!(" {name}=\"");
        let start = tag.find(&key).unwrap() + key.len();
        let len = tag[start..].find('"').unwrap();
        tag[start..start + len].parse().unwrap()
    }
    let tag_at = |i: usize| &doc[i..i + doc[i..].find('>').unwrap()];
    let count = doc
        .find(&format!("<{list} "))
        .map(|i| attr(tag_at(i), "count"));
    let indexes = doc
        .match_indices(&format!("<{elem} "))
        .map(|(i, _)| attr(tag_at(i), "index"))
        .collect();
    (count, indexes)
}

#[test]
fn document_with_skipped_records_validates_against_xsd() {
    if !xmllint_available() {
        eprintln!("skipping: xmllint not installed");
        return;
    }
    let opts = WriteOptions::new().with_record_error_policy(RecordErrorPolicy::Skip);
    for (more_errors, no_hint) in [(false, false), (true, false), (true, true)] {
        let mut src = Source::new("skip.raw", "xsd-test");
        src.with_error = true;
        src.with_more_errors = more_errors;
        src.no_hint = no_hint;
        // Source positions with the failing records counted, as a reader
        // would number them.
        src.spectra[1].index = 2;
        if more_errors {
            src.chroms[0].index = 1;
        }

        let mut plain = Vec::new();
        let summary = write_mzml_with_options(&mut src, &mut plain, &opts).unwrap();
        assert_eq!(summary.spectra_written, 2);
        assert_eq!(summary.spectra_skipped, 1 + usize::from(more_errors));
        assert_eq!(summary.chromatograms_skipped, usize::from(more_errors));
        validate(&plain, "mzML1.1.0.xsd", "xsd_skip_plain.mzML");

        let mut indexed = Vec::new();
        write_indexed_mzml_with_options(&mut src, &mut indexed, &opts).unwrap();
        validate(&indexed, "mzML1.1.2_idx.xsd", "xsd_skip_indexed.mzML");

        for doc in [&plain, &indexed] {
            let doc = std::str::from_utf8(doc).unwrap();
            assert_eq!(
                count_and_indexes(doc, "spectrumList", "spectrum"),
                (Some(2), vec![0, 1])
            );
            // The failing chromatogram (when present) comes first, so the
            // written one must be renumbered to 0.
            assert_eq!(
                count_and_indexes(doc, "chromatogramList", "chromatogram"),
                (Some(1), vec![0])
            );
        }
    }
}
