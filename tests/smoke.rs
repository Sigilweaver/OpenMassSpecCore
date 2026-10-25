//! End-to-end smoke test: build a minimal `SpectrumSource` the way a reader
//! crate would (constructors and setters only), write mzML, and check the
//! output contains the expected spectra.

use openmassspec_core::{
    write_indexed_mzml, write_mzml, Activation, Analyzer, CvTerm, MsPower, Polarity, PrecursorInfo,
    RunMetadata, ScanMode, SpectrumIter, SpectrumRecord, SpectrumSource,
};

struct ToySource {
    meta: RunMetadata,
    spectra: Vec<SpectrumRecord>,
    cursor: usize,
}

impl ToySource {
    fn new() -> Self {
        let meta = RunMetadata::new(
            "toy.raw",
            CvTerm::new("MS:1000563", "Thermo RAW format"),
            CvTerm::new("MS:1000768", "Thermo nativeID format"),
            CvTerm::new("MS:1001911", "Q Exactive"),
            "toy-writer",
            "0.0.0",
        )
        .with_analyzers(vec![Analyzer::FTMS]);
        let s1 = SpectrumRecord::new(0, 1, "controllerType=0 controllerNumber=1 scan=1")
            .with_ms_level(MsPower::Ms1.ms_level())
            .with_polarity(Polarity::Positive)
            .with_scan_mode(ScanMode::Centroid)
            .with_analyzer(Analyzer::FTMS)
            .with_filter("FTMS + p ESI Full ms".to_string())
            .with_retention_time_sec(0.123 * 60.0)
            .with_ion_injection_time_ms(20.0)
            .with_peaks(vec![100.0, 200.0, 300.0], vec![1.0, 5.0, 2.0]);
        let s2 = SpectrumRecord::new(1, 2, "controllerType=0 controllerNumber=1 scan=2")
            .with_ms_level(MsPower::Ms2.ms_level())
            .with_polarity(Polarity::Positive)
            .with_scan_mode(ScanMode::Centroid)
            .with_analyzer(Analyzer::FTMS)
            .with_filter("FTMS + p ESI d Full ms2 200.00@hcd28.00".to_string())
            .with_retention_time_sec(0.5 * 60.0)
            .with_total_ion_current(123.45)
            .with_base_peak_mz(150.5)
            .with_base_peak_intensity(99.0)
            .with_low_mz(100.0)
            .with_high_mz(180.0)
            .with_ion_injection_time_ms(50.0)
            .with_faims_cv(-45.0)
            .with_precursor(
                PrecursorInfo::default()
                    .with_target_mz(200.0)
                    .with_selected_mz(200.001)
                    .with_isolation_width(2.0)
                    .with_charge(2)
                    .with_collision_energy(28.0)
                    .with_ce_is_nce(true)
                    .with_precursor_native_id(
                        "controllerType=0 controllerNumber=1 scan=1".to_string(),
                    )
                    .with_activation(Activation::CID)
                    .with_analyzer(Analyzer::FTMS),
            )
            .with_peaks(vec![150.5, 160.0], vec![99.0, 50.0]);
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

    fn iter_spectra<'a>(&'a mut self) -> SpectrumIter<'a> {
        self.cursor = 0;
        Box::new(std::iter::from_fn(move || {
            let rec = self.spectra.get(self.cursor)?.clone();
            self.cursor += 1;
            Some(Ok(rec))
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
    let summary = write_mzml(&mut src, &mut buf).unwrap();
    assert_eq!(summary.spectra_written, 2);
    let s = String::from_utf8(buf).unwrap();

    assert!(s.starts_with(r#"<?xml version="1.0" encoding="utf-8"?>"#));
    assert!(s.contains(r#"<mzML xmlns="http://psi.hupo.org/ms/mzml""#));
    assert!(s.contains(r#"<spectrumList count="2""#));
    assert!(s.contains(r#"id="controllerType=0 controllerNumber=1 scan=1""#));
    assert!(s.contains(r#"id="controllerType=0 controllerNumber=1 scan=2""#));
    assert!(s.contains(r#"<cvParam cvRef="MS" accession="MS:1001911" name="Q Exactive""#));
    // CID stays resonant CID whatever analyzer read the fragments out.
    assert!(s.contains(r#"accession="MS:1000133" name="collision-induced dissociation""#));
    assert!(!s.contains("MS:1000422"));
    // s2 carries a FAIMS CV; s1 does not.
    assert!(
        s.contains(r#"accession="MS:1001581" name="FAIMS compensation voltage" value="-45.0000""#)
    );
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
    let recs: Vec<_> = src
        .iter_spectra()
        .collect::<openmassspec_core::Result<_>>()
        .unwrap();
    assert_eq!(recs[0].effective_tic(), 8.0);
    let (bp_mz, bp_i) = recs[0].effective_base_peak().unwrap();
    assert!((bp_mz - 200.0).abs() < 1e-9);
    assert!((bp_i - 5.0).abs() < 1e-9);
    let (lo, hi) = recs[0].effective_mz_range().unwrap();
    assert_eq!(lo, 100.0);
    assert_eq!(hi, 300.0);
}
