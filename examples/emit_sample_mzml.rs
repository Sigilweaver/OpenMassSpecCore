//! Emit a small synthetic mzML document (plain and indexed) from the
//! canonical writer, for schema-conformance checking in CI.
//!
//! Usage:
//!   emit_sample_mzml [PLAIN_OUT] [INDEXED_OUT]
//!
//! Defaults to `sample_plain.mzML` and `sample_indexed.mzML` in the
//! current directory. The output is validated against the vendored
//! PSI-MS mzML XSD by the `validate-mzml` CI job.

use std::fs::File;
use std::io::BufWriter;

use openmassspec_core::{
    write_indexed_mzml, write_mzml, Activation, Analyzer, CvTerm, MobilityArrayKind, MsPower,
    Polarity, PrecursorInfo, RunMetadata, ScanMode, SpectrumRecord, SpectrumSource,
};

struct SampleSource {
    meta: RunMetadata,
    spectra: Vec<SpectrumRecord>,
    cursor: usize,
}

impl SampleSource {
    fn new() -> Self {
        let meta = RunMetadata {
            source_file_name: "sample.raw".into(),
            source_file_format: CvTerm::new("MS:1000563", "Thermo RAW format"),
            native_id_format: CvTerm::new("MS:1000768", "Thermo nativeID format"),
            instrument: CvTerm::new("MS:1001911", "Q Exactive"),
            software_name: "openmassspec-core-sample".into(),
            software_version: env!("CARGO_PKG_VERSION").into(),
            start_timestamp: None,
            mobility_array_kind: Some(MobilityArrayKind::InverseReducedVsPerCm2),
        };
        let ms1 = SpectrumRecord {
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
        let ms2 = SpectrumRecord {
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
        // A frame-collapsed ion-mobility MS1 (the Bruker timsTOF shape):
        // scalar 1/K0 plus a per-peak inverse-reduced-mobility array, which
        // the writer emits as a third binary data array. Exercises the
        // ion-mobility CV path so it is schema-checked without a vendor file.
        let ms1_mobility = SpectrumRecord {
            index: 2,
            scan_number: 3,
            native_id: "frame=1 scan=0".into(),
            ms_level: MsPower::Ms1.ms_level(),
            polarity: Some(Polarity::Positive),
            scan_mode: Some(ScanMode::Centroid),
            analyzer: Some(Analyzer::TOFMS),
            filter: None,
            retention_time_sec: 0.75 * 60.0,
            total_ion_current: None,
            base_peak_mz: None,
            base_peak_intensity: None,
            low_mz: None,
            high_mz: None,
            ion_injection_time_ms: None,
            inv_mobility: Some(0.95),
            precursor: None,
            mz: vec![120.0, 240.0, 360.0],
            intensity: vec![3.0, 7.0, 4.0],
            inv_mobility_per_peak: Some(vec![0.92, 0.95, 0.98]),
        };
        Self {
            meta,
            spectra: vec![ms1, ms2, ms1_mobility],
            cursor: 0,
        }
    }
}

impl SpectrumSource for SampleSource {
    fn run_metadata(&self) -> RunMetadata {
        self.meta.clone()
    }

    fn iter_spectra<'a>(&'a mut self) -> Box<dyn Iterator<Item = SpectrumRecord> + 'a> {
        self.cursor = 0;
        Box::new(std::iter::from_fn(move || {
            let rec = self.spectra.get(self.cursor).cloned();
            if rec.is_some() {
                self.cursor += 1;
            }
            rec
        }))
    }

    fn spectrum_count_hint(&self) -> Option<usize> {
        Some(self.spectra.len())
    }
}

fn main() -> std::io::Result<()> {
    let mut args = std::env::args().skip(1);
    let plain_path = args.next().unwrap_or_else(|| "sample_plain.mzML".into());
    let indexed_path = args.next().unwrap_or_else(|| "sample_indexed.mzML".into());

    let mut src = SampleSource::new();
    let mut plain = BufWriter::new(File::create(&plain_path)?);
    write_mzml(&mut src, &mut plain)
        .map_err(|e| std::io::Error::other(format!("write_mzml: {e}")))?;

    let mut indexed = BufWriter::new(File::create(&indexed_path)?);
    write_indexed_mzml(&mut src, &mut indexed)
        .map_err(|e| std::io::Error::other(format!("write_indexed_mzml: {e}")))?;

    eprintln!("wrote {plain_path} and {indexed_path}");
    Ok(())
}
