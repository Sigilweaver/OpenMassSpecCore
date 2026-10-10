//! Validates writer output against the vendored PSI-MS mzML XSDs with
//! `xmllint`, covering shapes the CI `validate-mzml` example does not: a
//! digit-leading source file name (the `<run id>` is an `xs:ID`), an MS2
//! precursor with unknown activation, and an SRM chromatogram precursor
//! (both of which now write an empty `<activation>`).
//!
//! Skips with a message when `xmllint` is not installed.

use std::path::{Path, PathBuf};
use std::process::Command;

use openmassspec_core::{
    write_indexed_mzml, write_mzml, ChromatogramRecord, CvTerm, Polarity, PrecursorInfo,
    RunMetadata, ScanMode, SpectrumRecord, SpectrumSource,
};

struct Source {
    meta: RunMetadata,
    spectra: Vec<SpectrumRecord>,
    chroms: Vec<ChromatogramRecord>,
}

fn spectrum(index: usize, ms_level: u32, precursor: Option<PrecursorInfo>) -> SpectrumRecord {
    SpectrumRecord {
        index,
        scan_number: index as u32 + 1,
        native_id: format!("controllerType=0 controllerNumber=1 scan={}", index + 1),
        ms_level,
        polarity: Some(Polarity::Positive),
        scan_mode: Some(ScanMode::Centroid),
        analyzer: None,
        acquisition_event_id: None,
        filter: None,
        retention_time_sec: index as f64,
        total_ion_current: None,
        base_peak_mz: None,
        base_peak_intensity: None,
        low_mz: None,
        high_mz: None,
        ion_injection_time_ms: None,
        inv_mobility: None,
        faims_cv: None,
        precursor,
        mz: vec![100.0, 200.0],
        intensity: vec![1.0, 2.0],
        inv_mobility_per_peak: None,
        extra: Default::default(),
    }
}

impl Source {
    fn new(source_file_name: &str) -> Self {
        Self {
            meta: RunMetadata {
                source_file_name: source_file_name.into(),
                source_file_format: CvTerm::new("MS:1000563", "Thermo RAW format"),
                native_id_format: CvTerm::new("MS:1000768", "Thermo nativeID format"),
                instrument: CvTerm::new("MS:1001911", "Q Exactive"),
                instrument_serial_number: None,
                acquisition_software_name: Some("Xcalibur".into()),
                acquisition_software_version: Some("4.4".into()),
                software_name: "xsd-test".into(),
                software_version: "0.0.0".into(),
                start_timestamp: None,
                mobility_array_kind: None,
                analyzers: Vec::new(),
                extra: Default::default(),
            },
            spectra: vec![
                spectrum(0, 1, None),
                spectrum(
                    1,
                    2,
                    Some(PrecursorInfo {
                        selected_mz: Some(150.0),
                        collision_energy: Some(30.0),
                        ce_is_nce: true,
                        activation: None,
                        ..Default::default()
                    }),
                ),
            ],
            chroms: vec![ChromatogramRecord {
                index: 0,
                id: "SRM SIC Q1=524.3 Q3=136.1".into(),
                chromatogram_type: Some(CvTerm::new(
                    "MS:1001473",
                    "selected reaction monitoring chromatogram",
                )),
                precursor_mz: Some(524.3),
                product_mz: Some(136.1),
                time_sec: vec![0.0, 30.0],
                intensity: vec![5.0, 8.0],
            }],
        }
    }
}

impl SpectrumSource for Source {
    fn run_metadata(&self) -> RunMetadata {
        self.meta.clone()
    }
    fn iter_spectra<'a>(&'a mut self) -> Box<dyn Iterator<Item = SpectrumRecord> + 'a> {
        Box::new(self.spectra.clone().into_iter())
    }
    fn iter_chromatograms<'a>(&'a mut self) -> Box<dyn Iterator<Item = ChromatogramRecord> + 'a> {
        Box::new(self.chroms.clone().into_iter())
    }
    fn spectrum_count_hint(&self) -> Option<usize> {
        Some(self.spectra.len())
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
fn digit_leading_file_name_validates_against_xsd() {
    if !xmllint_available() {
        eprintln!("skipping: xmllint not installed");
        return;
    }
    for name in ["20260101_sample.raw", "1 run (copy).raw"] {
        let mut src = Source::new(name);

        let mut plain = Vec::new();
        write_mzml(&mut src, &mut plain).unwrap();
        validate(&plain, "mzML1.1.0.xsd", "xsd_plain.mzML");

        let mut indexed = Vec::new();
        write_indexed_mzml(&mut src, &mut indexed).unwrap();
        validate(&indexed, "mzML1.1.2_idx.xsd", "xsd_indexed.mzML");
    }
}
