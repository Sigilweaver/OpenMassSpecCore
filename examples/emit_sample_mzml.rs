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
use std::io::{BufWriter, Write};

use openmassspec_core::{
    write_indexed_mzml, write_mzml, Activation, Analyzer, ChromatogramIter, ChromatogramRecord,
    CvTerm, MobilityArrayKind, MsPower, Polarity, PrecursorInfo, RunMetadata, ScanMode,
    SpectrumIter, SpectrumRecord, SpectrumSource,
};

struct SampleSource {
    meta: RunMetadata,
    spectra: Vec<SpectrumRecord>,
    chroms: Vec<ChromatogramRecord>,
}

impl SampleSource {
    fn new() -> Self {
        let meta = RunMetadata::new(
            // Digit-leading on purpose: `<run id>` is an `xs:ID`, so this
            // checks the writer turns it into a valid NCName.
            "20260601_sample.raw",
            CvTerm::new("MS:1000563", "Thermo RAW format"),
            CvTerm::new("MS:1000768", "Thermo nativeID format"),
            CvTerm::new("MS:1001911", "Q Exactive"),
            // Not an NCName on purpose: `<software id>` is an `xs:ID` too.
            "openmassspec-core sample",
            env!("CARGO_PKG_VERSION"),
        )
        .with_instrument_serial_number("SN-EXAMPLE-001".to_string())
        // Exercises the second `<software>` entry / `<softwareList
        // count="2">` path so it is schema-checked without a vendor file.
        .with_acquisition_software_name("Xcalibur".to_string())
        .with_acquisition_software_version("4.4.16.14".to_string())
        // Exercises the `startTimeStamp` attribute path.
        .with_start_timestamp("2026-06-01T14:30:00Z".to_string())
        .with_mobility_array_kind(MobilityArrayKind::InverseReducedVsPerCm2)
        // Exercises the per-spectrum `instrumentConfigurationRef` path:
        // the FTMS (Orbitrap) spectra below get one instrumentConfiguration
        // and the TOFMS one gets another.
        .with_analyzers(vec![Analyzer::FTMS, Analyzer::TOFMS])
        .with_extra("openmassspec.sample_source", "synthetic");

        let ms1 = SpectrumRecord::new(0, 1, "controllerType=0 controllerNumber=1 scan=1")
            .with_ms_level(MsPower::Ms1.ms_level())
            .with_polarity(Polarity::Positive)
            .with_scan_mode(ScanMode::Centroid)
            .with_analyzer(Analyzer::FTMS)
            .with_acquisition_event_id(1)
            .with_filter("FTMS + p ESI Full ms".to_string())
            .with_retention_time_sec(0.123 * 60.0)
            .with_ion_injection_time_ms(20.0)
            .with_peaks(vec![100.0, 200.0, 300.0], vec![1.0, 5.0, 2.0])
            .with_extra("openmassspec.sample_label", "survey");
        let ms2 = SpectrumRecord::new(1, 2, "controllerType=0 controllerNumber=1 scan=2")
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
                    .with_collision_energy(28.0, openmassspec_core::CollisionEnergyUnit::Nce)
                    .with_precursor_native_id(
                        "controllerType=0 controllerNumber=1 scan=1".to_string(),
                    )
                    .with_activation(Activation::CID)
                    .with_analyzer(Analyzer::FTMS)
                    // Exercises the `collisional cross sectional area`
                    // selectedIon cvParam path.
                    .with_ccs(153.4),
            )
            .with_peaks(vec![150.5, 160.0], vec![99.0, 50.0]);
        // A frame-collapsed ion-mobility MS1 (the Bruker timsTOF shape):
        // scalar 1/K0 plus a per-peak inverse-reduced-mobility array, which
        // the writer emits as a third binary data array.
        let ms1_mobility = SpectrumRecord::new(2, 3, "frame=1 scan=0")
            .with_ms_level(MsPower::Ms1.ms_level())
            .with_polarity(Polarity::Positive)
            .with_scan_mode(ScanMode::Centroid)
            .with_analyzer(Analyzer::TOFMS)
            .with_retention_time_sec(0.75 * 60.0)
            .with_inv_mobility(0.95)
            .with_peaks(vec![120.0, 240.0, 360.0], vec![3.0, 7.0, 4.0])
            .with_inv_mobility_per_peak(vec![0.92, 0.95, 0.98]);
        // A TIC and a SRM/MRM transition chromatogram. Exercises the
        // `chromatogramList` path (plain + indexed).
        let tic =
            ChromatogramRecord::new(0, "TIC", vec![0.0, 30.0, 45.0], vec![120.0, 340.0, 210.0])
                .with_chromatogram_type(CvTerm::new(
                    "MS:1000235",
                    "total ion current chromatogram",
                ));
        let srm = ChromatogramRecord::new(
            1,
            "SRM SIC Q1=524.3 Q3=136.1",
            vec![0.0, 30.0],
            vec![50.0, 80.0],
        )
        .with_chromatogram_type(CvTerm::new(
            "MS:1001473",
            "selected reaction monitoring chromatogram",
        ))
        .with_precursor_mz(524.3)
        .with_product_mz(136.1);
        Self {
            meta,
            spectra: vec![ms1, ms2, ms1_mobility],
            chroms: vec![tic, srm],
        }
    }
}

impl SpectrumSource for SampleSource {
    fn run_metadata(&self) -> RunMetadata {
        self.meta.clone()
    }

    fn iter_spectra<'a>(&'a mut self) -> SpectrumIter<'a> {
        Box::new(self.spectra.iter().cloned().map(Ok))
    }

    fn iter_chromatograms<'a>(&'a mut self) -> ChromatogramIter<'a> {
        Box::new(self.chroms.iter().cloned().map(Ok))
    }

    fn spectrum_count_hint(&self) -> Option<usize> {
        Some(self.spectra.len())
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let plain_path = args.next().unwrap_or_else(|| "sample_plain.mzML".into());
    let indexed_path = args.next().unwrap_or_else(|| "sample_indexed.mzML".into());

    let mut src = SampleSource::new();
    let mut plain = BufWriter::new(File::create(&plain_path)?);
    write_mzml(&mut src, &mut plain)?;
    plain.flush()?;

    let mut indexed = BufWriter::new(File::create(&indexed_path)?);
    write_indexed_mzml(&mut src, &mut indexed)?;
    indexed.flush()?;

    eprintln!("wrote {plain_path} and {indexed_path}");
    Ok(())
}
