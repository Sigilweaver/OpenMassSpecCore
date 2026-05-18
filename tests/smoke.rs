//! End-to-end smoke test: build a minimal `SpectrumSource`, write mzML,
//! check the output is well-formed XML and contains the expected spectra.

use mass_spec_core::{
    write_indexed_mzml, write_mzml, Activation, Analyzer, CvTerm, MsPower, Polarity, PrecursorInfo,
    RunMetadata, ScanMode, SpectrumRecord, SpectrumSource,
};

struct ToySource {
    meta: RunMetadata,
    spectra: Vec<SpectrumRecord>,
    cursor: usize,
}

impl ToySource {
    fn new() -> Self {
        let meta = RunMetadata {
            source_file_name: "toy.raw".into(),
            source_file_format: CvTerm::new("MS:1000563", "Thermo RAW format"),
            native_id_format: CvTerm::new("MS:1000768", "Thermo nativeID format"),
            instrument: CvTerm::new("MS:1001911", "Q Exactive"),
            software_name: "toy-writer".into(),
            software_version: "0.0.0".into(),
            start_timestamp: None,
        };
        let s1 = SpectrumRecord {
            index: 0,
            scan_number: 1,
            native_id: "controllerType=0 controllerNumber=1 scan=1".into(),
            ms_level: MsPower::Ms1.ms_level(),
            polarity: Some(Polarity::Positive),
            scan_mode: Some(ScanMode::Centroid),
            analyzer: Some(Analyzer::FTMS),
            filter: Some("FTMS + p ESI Full ms".into()),
            retention_time_sec: 0.123 * 60.0,
            total_ion_current: None,
            base_peak_mz: None,
            base_peak_intensity: None,
            low_mz: None,
            high_mz: None,
            ion_injection_time_ms: Some(20.0),
            inv_mobility: None,
            precursor: None,
            mz: vec![100.0, 200.0, 300.0],
            intensity: vec![1.0, 5.0, 2.0],
            inv_mobility_per_peak: None,
        };
        let s2 = SpectrumRecord {
            index: 1,
            scan_number: 2,
            native_id: "controllerType=0 controllerNumber=1 scan=2".into(),
            ms_level: MsPower::Ms2.ms_level(),
            polarity: Some(Polarity::Positive),
            scan_mode: Some(ScanMode::Centroid),
            analyzer: Some(Analyzer::FTMS),
            filter: Some("FTMS + p ESI d Full ms2 200.00@hcd28.00".into()),
            retention_time_sec: 0.5 * 60.0,
            total_ion_current: Some(123.45),
            base_peak_mz: Some(150.5),
            base_peak_intensity: Some(99.0),
            low_mz: Some(100.0),
            high_mz: Some(180.0),
            ion_injection_time_ms: Some(50.0),
            inv_mobility: None,
            precursor: Some(PrecursorInfo {
                target_mz: Some(200.0),
                selected_mz: Some(200.001),
                isolation_width: Some(2.0),
                charge: Some(2),
                intensity: None,
                collision_energy: Some(28.0),
                ce_is_nce: true,
                precursor_native_id: Some("controllerType=0 controllerNumber=1 scan=1".into()),
                activation: Some(Activation::CID),
                analyzer: Some(Analyzer::FTMS),
            }),
            mz: vec![150.5, 160.0],
            intensity: vec![99.0, 50.0],
            inv_mobility_per_peak: None,
        };
        Self {
            meta,
            spectra: vec![s1, s2],
            cursor: 0,
        }
    }
}

impl SpectrumSource for ToySource {
    fn run_metadata(&self) -> RunMetadata {
        self.meta.clone()
    }

    fn iter_spectra<'a>(&'a mut self) -> Box<dyn Iterator<Item = SpectrumRecord> + 'a> {
        self.cursor = 0;
        Box::new(std::iter::from_fn(move || {
            if self.cursor >= self.spectra.len() {
                None
            } else {
                let rec = self.spectra[self.cursor].clone();
                self.cursor += 1;
                Some(rec)
            }
        }))
    }

    fn spectrum_count_hint(&self) -> Option<usize> {
        Some(self.spectra.len())
    }
}

#[test]
fn writes_plain_mzml() {
    let mut src = ToySource::new();
    let mut buf = Vec::new();
    write_mzml(&mut src, &mut buf).unwrap();
    let s = String::from_utf8(buf).unwrap();

    assert!(s.starts_with(r#"<?xml version="1.0" encoding="utf-8"?>"#));
    assert!(s.contains(r#"<mzML xmlns="http://psi.hupo.org/ms/mzml""#));
    assert!(s.contains(r#"<spectrumList count="2""#));
    assert!(s.contains(r#"id="controllerType=0 controllerNumber=1 scan=1""#));
    assert!(s.contains(r#"id="controllerType=0 controllerNumber=1 scan=2""#));
    assert!(s.contains(r#"<cvParam cvRef="MS" accession="MS:1001911" name="Q Exactive""#));
    // CID on FTMS analyzer should map to beam-type CID, not ion-trap CID.
    assert!(s.contains(r#"accession="MS:1000422" name="beam-type collision-induced dissociation""#));
    assert!(s.ends_with("</mzML>\n"));
}

#[test]
fn writes_indexed_mzml() {
    let mut src = ToySource::new();
    let mut buf = Vec::new();
    write_indexed_mzml(&mut src, &mut buf).unwrap();
    let s = String::from_utf8(buf).unwrap();

    assert!(s.starts_with(r#"<?xml version="1.0" encoding="utf-8"?>"#));
    assert!(s.contains("<indexedmzML"));
    assert!(s.contains("<indexList count=\"1\""));
    assert!(s.contains("<indexListOffset>"));
    assert!(s.contains("<fileChecksum>"));
    assert!(s.contains(r#"idRef="controllerType=0 controllerNumber=1 scan=1""#));
    assert!(s.ends_with("</indexedmzML>\n"));
}

#[test]
fn effective_helpers_compute_from_arrays() {
    let mut src = ToySource::new();
    // s1 has no pre-populated TIC; should sum to 8.0.
    let recs: Vec<_> = src.iter_spectra().collect();
    assert_eq!(recs[0].effective_tic(), 8.0);
    let (bp_mz, bp_i) = recs[0].effective_base_peak().unwrap();
    assert!((bp_mz - 200.0).abs() < 1e-9);
    assert!((bp_i - 5.0).abs() < 1e-9);
    let (lo, hi) = recs[0].effective_mz_range().unwrap();
    assert_eq!(lo, 100.0);
    assert_eq!(hi, 300.0);
}
