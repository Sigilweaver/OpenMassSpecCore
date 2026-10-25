//! Vendor-neutral profile-to-centroid peak picking.
//!
//! Centroiding is a transform over the arrays a [`SpectrumSource`] already
//! yields (`mz`/`intensity`, and optionally `inv_mobility_per_peak`) - once a
//! [`SpectrumRecord`] exists, the operation is identical regardless of which
//! vendor parser produced it. [`Centroided`] wraps any `SpectrumSource` and
//! applies the transform lazily, one spectrum at a time, so it composes with
//! the streaming mzML writer and the Arrow bridge without buffering a whole
//! run in memory.

use crate::cv;
use crate::enums::ScanMode;
use crate::error::{Error, RecordProblem, Result};
use crate::leniency::{Check, Leniency};
use crate::source::{ChromatogramIter, SpectrumIter, SpectrumSource};
use crate::types::{CvTerm, RunMetadata, SpectrumRecord};
use crate::validate;
use crate::warning::{Warning, WarningKind};

/// Relative median m/z spacing above which a spectrum labelled profile is
/// taken to be a peak list (0.5 percent of m/z, i.e. 5000 ppm). Profile
/// sampling is far finer than that on every analyzer family (Orbitrap and
/// TOF are tens of ppm; low-resolution ion-trap steps are hundreds).
const PEAK_LIST_RELATIVE_SPACING: f64 = 0.005;

/// A neighbour joins a peak only if it is within this many median spacings.
const NEIGHBOR_SPACING_FACTOR: f64 = 3.0;

/// Wraps a [`SpectrumSource`], centroiding every spectrum declared
/// [`ScanMode::Profile`]. Spectra declared [`ScanMode::Centroid`] pass
/// through unchanged (idempotent).
///
/// A spectrum with no declared scan mode is a record error
/// ([`RecordProblem::UnknownScanMode`]) unless
/// [`with_leniency`](Centroided::with_leniency) allows
/// [`Check::UnknownScanMode`], which passes it through untouched with a
/// [`WarningKind::CentroidNotApplied`] warning. A spectrum declared
/// profile whose points are too sparse to be profile data is a record error
/// ([`RecordProblem::NotProfile`]). Peak weighting only joins neighbours
/// that sit within a few median sample spacings, so gaps in a
/// zero-stripped profile are never averaged across.
pub struct Centroided<S: SpectrumSource> {
    inner: S,
    min_intensity: f32,
    leniency: Leniency,
}

impl<S: SpectrumSource> Centroided<S> {
    /// Wrap `inner`, picking every local-maximum peak regardless of height.
    pub fn new(inner: S) -> Self {
        Self {
            inner,
            min_intensity: 0.0,
            leniency: Leniency::strict(),
        }
    }

    /// Discard picked peaks below `min_intensity` (a simple noise floor).
    pub fn with_min_intensity(mut self, min_intensity: f32) -> Self {
        self.min_intensity = min_intensity;
        self
    }

    /// Checks that may continue with a warning instead of failing.
    pub fn with_leniency(mut self, leniency: Leniency) -> Self {
        self.leniency = leniency;
        self
    }
}

impl<S: SpectrumSource> SpectrumSource for Centroided<S> {
    fn run_metadata(&self) -> RunMetadata {
        self.inner.run_metadata()
    }

    /// Errors from the wrapped source pass through unchanged.
    fn iter_spectra<'a>(&'a mut self) -> SpectrumIter<'a> {
        let min_intensity = self.min_intensity;
        let leniency = self.leniency.clone();
        Box::new(
            self.inner
                .iter_spectra()
                .map(move |rec| rec.and_then(|rec| centroid_record(rec, min_intensity, &leniency))),
        )
    }

    fn iter_chromatograms<'a>(&'a mut self) -> ChromatogramIter<'a> {
        self.inner.iter_chromatograms()
    }

    fn spectrum_count_hint(&self) -> Option<usize> {
        self.inner.spectrum_count_hint()
    }

    fn additional_processing_steps(&self) -> Vec<CvTerm> {
        let mut steps = self.inner.additional_processing_steps();
        steps.push(cv::PEAK_PICKING.into());
        steps
    }
}

fn record_error(rec: &SpectrumRecord, problem: RecordProblem) -> Error {
    Error::Record {
        record_id: rec.native_id.clone(),
        problem,
    }
}

fn centroid_record(
    mut rec: SpectrumRecord,
    min_intensity: f32,
    leniency: &Leniency,
) -> Result<SpectrumRecord> {
    match rec.scan_mode {
        Some(ScanMode::Centroid) => return Ok(rec),
        Some(ScanMode::Profile) => {}
        _ => {
            if !leniency.allows(Check::UnknownScanMode) {
                return Err(record_error(&rec, RecordProblem::UnknownScanMode));
            }
            rec.warnings.push(Warning::new(
                WarningKind::CentroidNotApplied,
                "scan mode unknown: left as read, not centroided",
            ));
            return Ok(rec);
        }
    }

    validate::check_array_lengths(&rec).map_err(|p| record_error(&rec, p))?;
    if let Some((spacing, relative)) = validate::median_spacing(&rec.mz) {
        if relative > PEAK_LIST_RELATIVE_SPACING {
            return Err(record_error(
                &rec,
                RecordProblem::NotProfile {
                    median_spacing: spacing,
                    relative_ppm: relative * 1e6,
                },
            ));
        }
    }

    let (mz, intensity, inv_mobility_per_peak) = pick_peaks(
        &rec.mz,
        &rec.intensity,
        rec.inv_mobility_per_peak.as_deref(),
        min_intensity,
    );

    rec.mz = mz;
    rec.intensity = intensity;
    rec.inv_mobility_per_peak = inv_mobility_per_peak;
    rec.scan_mode = Some(ScanMode::Centroid);
    // These were derived from the profile arrays; let the mzML writer /
    // Arrow bridge recompute them from the new centroided arrays instead of
    // carrying stale values forward (see `SpectrumRecord`'s field docs).
    rec.total_ion_current = None;
    rec.base_peak_mz = None;
    rec.base_peak_intensity = None;
    rec.low_mz = None;
    rec.high_mz = None;
    Ok(rec)
}

/// Local-maxima peak picking over contiguous runs of samples. A neighbour
/// counts only if it is within `NEIGHBOR_SPACING_FACTOR` median spacings; a
/// farther one is a gap in a zero-stripped profile and is treated as absent.
/// A point is a picked peak if it is no smaller than both counted neighbors
/// and strictly larger than at least one of them (a point with no counted
/// neighbor is a peak on its own). The picked m/z (and, when present,
/// inverse mobility) is the intensity-weighted centroid over the apex and
/// its counted neighbors; the picked intensity is the apex height. This is
/// intentionally simple - see `timsrust-centroid` / `pyteomics` for more
/// sophisticated peer approaches if this ever needs to improve on plain
/// local-maxima picking.
fn pick_peaks(
    mz: &[f64],
    intensity: &[f32],
    inv_mobility_per_peak: Option<&[f32]>,
    min_intensity: f32,
) -> (Vec<f64>, Vec<f32>, Option<Vec<f32>>) {
    let n = mz.len();
    let mut out_mz = Vec::new();
    let mut out_intensity = Vec::new();
    let mut out_im = inv_mobility_per_peak.map(|_| Vec::new());

    let gap_bound = match validate::median_spacing(mz) {
        Some((spacing, _)) if spacing > 0.0 => spacing * NEIGHBOR_SPACING_FACTOR,
        _ => {
            let mid = mz.get(n / 2).copied().unwrap_or(0.0);
            mid.abs() * PEAK_LIST_RELATIVE_SPACING
        }
    };
    let has_left = |i: usize| i > 0 && mz[i] - mz[i - 1] <= gap_bound;
    let has_right = |i: usize| i + 1 < n && mz[i + 1] - mz[i] <= gap_bound;

    for i in 0..n {
        let left = has_left(i);
        let right = has_right(i);
        let no_smaller_than_left = !left || intensity[i] >= intensity[i - 1];
        let no_smaller_than_right = !right || intensity[i] >= intensity[i + 1];
        let beats_a_neighbor = (!left && !right)
            || (left && intensity[i] > intensity[i - 1])
            || (right && intensity[i] > intensity[i + 1]);
        if !(no_smaller_than_left && no_smaller_than_right && beats_a_neighbor) {
            continue;
        }
        if intensity[i] < min_intensity {
            continue;
        }

        let lo = if left { i - 1 } else { i };
        let hi = if right { i + 1 } else { i };
        let mut weighted_mz = 0.0f64;
        let mut weighted_im = 0.0f64;
        let mut weight_sum = 0.0f64;
        for j in lo..=hi {
            let w = f64::from(intensity[j]);
            weighted_mz += mz[j] * w;
            weight_sum += w;
            if let Some(im) = inv_mobility_per_peak {
                weighted_im += f64::from(im[j]) * w;
            }
        }
        let centroid_mz = if weight_sum > 0.0 {
            weighted_mz / weight_sum
        } else {
            mz[i]
        };
        out_mz.push(centroid_mz);
        out_intensity.push(intensity[i]);
        if let Some(out) = out_im.as_mut() {
            #[allow(clippy::cast_possible_truncation)]
            let centroid_im = if weight_sum > 0.0 {
                (weighted_im / weight_sum) as f32
            } else {
                inv_mobility_per_peak.expect("out_im is Some only when inv_mobility_per_peak is")[i]
            };
            out.push(centroid_im);
        }
    }

    (out_mz, out_intensity, out_im)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformance::assert_source_invariants;
    use crate::{MsPower, Polarity};

    struct OneSpectrumSource {
        meta: RunMetadata,
        rec: Option<SpectrumRecord>,
        /// Yield an `Err` after the record.
        then_fail: bool,
    }

    impl OneSpectrumSource {
        fn new(rec: SpectrumRecord) -> Self {
            Self {
                meta: RunMetadata {
                    source_file_name: "toy.raw".into(),
                    source_file_format: CvTerm::new("MS:1000563", "Thermo RAW format"),
                    native_id_format: CvTerm::new("MS:1000768", "Thermo nativeID format"),
                    instrument: CvTerm::new("MS:1001911", "Q Exactive"),
                    instrument_serial_number: None,
                    software_name: "toy-writer".into(),
                    software_version: "0.0.0".into(),
                    acquisition_software_name: None,
                    acquisition_software_version: None,
                    start_timestamp: None,
                    mobility_array_kind: None,
                    analyzers: Vec::new(),
                    warnings: Vec::new(),
                    extra: Default::default(),
                },
                rec: Some(rec),
                then_fail: false,
            }
        }
    }

    impl SpectrumSource for OneSpectrumSource {
        fn run_metadata(&self) -> RunMetadata {
            self.meta.clone()
        }

        fn iter_spectra<'a>(&'a mut self) -> SpectrumIter<'a> {
            let fail = self.then_fail.then(|| {
                Err(crate::Error::vendor(std::io::Error::other(
                    "truncated scan",
                )))
            });
            Box::new(self.rec.take().map(Ok).into_iter().chain(fail))
        }

        fn spectrum_count_hint(&self) -> Option<usize> {
            Some(usize::from(self.rec.is_some()))
        }
    }

    fn profile_spectrum() -> SpectrumRecord {
        // Two synthetic Gaussian-ish bumps centered at mz=100 and mz=200,
        // known peak apexes so we can assert the picked centroids land near
        // them. Compared only against this hand-built synthetic input, not
        // any vendor tool or reference converter (clean-room rule).
        let mz: Vec<f64> = vec![
            99.0, 99.5, 100.0, 100.5, 101.0, // bump 1, apex at 100.0
            199.0, 199.5, 200.0, 200.5, 201.0, // bump 2, apex at 200.0
        ];
        let intensity: Vec<f32> = vec![10.0, 50.0, 100.0, 50.0, 10.0, 5.0, 40.0, 80.0, 40.0, 5.0];
        SpectrumRecord {
            index: 0,
            scan_number: 1,
            native_id: "scan=1".into(),
            ms_level: MsPower::Ms1.ms_level(),
            polarity: Some(Polarity::Positive),
            scan_mode: Some(ScanMode::Profile),
            analyzer: None,
            acquisition_event_id: None,
            filter: None,
            retention_time_sec: 1.0,
            total_ion_current: Some(intensity.iter().map(|&v| v as f64).sum()),
            base_peak_mz: Some(100.0),
            base_peak_intensity: Some(100.0),
            low_mz: Some(99.0),
            high_mz: Some(201.0),
            scan_window: None,
            ion_injection_time_ms: None,
            inv_mobility: None,
            faims_cv: None,
            precursor: None,
            mz,
            intensity,
            inv_mobility_per_peak: None,
            warnings: Vec::new(),
            extra: Default::default(),
        }
    }

    #[test]
    fn centroids_profile_spectrum_at_known_apexes() {
        let mut src = Centroided::new(OneSpectrumSource::new(profile_spectrum()));
        let recs: Vec<_> = src.iter_spectra().collect::<crate::Result<_>>().unwrap();
        assert_eq!(recs.len(), 1);
        let rec = &recs[0];

        assert_eq!(rec.scan_mode, Some(ScanMode::Centroid));
        assert!(
            rec.mz.len() < 10,
            "expected fewer peaks after centroiding, got {}",
            rec.mz.len()
        );
        assert_eq!(rec.mz.len(), rec.intensity.len());

        let near_100 = rec.mz.iter().any(|&m| (m - 100.0).abs() < 0.6);
        let near_200 = rec.mz.iter().any(|&m| (m - 200.0).abs() < 0.6);
        assert!(near_100, "no picked peak near mz=100: {:?}", rec.mz);
        assert!(near_200, "no picked peak near mz=200: {:?}", rec.mz);

        // Stale profile-derived summary fields must be cleared, not carried
        // forward, so the writer recomputes them from the new arrays.
        assert_eq!(rec.total_ion_current, None);
        assert_eq!(rec.base_peak_mz, None);
    }

    #[test]
    fn already_centroided_spectrum_passes_through_unchanged() {
        let mut rec = profile_spectrum();
        rec.scan_mode = Some(ScanMode::Centroid);
        let original_mz = rec.mz.clone();

        let mut src = Centroided::new(OneSpectrumSource::new(rec));
        let recs: Vec<_> = src.iter_spectra().collect::<crate::Result<_>>().unwrap();
        assert_eq!(recs[0].mz, original_mz);
    }

    #[test]
    fn additional_processing_steps_reports_peak_picking() {
        let src = Centroided::new(OneSpectrumSource::new(profile_spectrum()));
        assert_eq!(
            src.additional_processing_steps(),
            vec![CvTerm::from(cv::PEAK_PICKING)]
        );
    }

    #[test]
    fn additional_processing_steps_accumulate_through_nested_adapters() {
        // Wrapping a Centroided source in another Centroided (a no-op in
        // practice, since it is idempotent) must not drop the inner
        // adapter's step - steps accumulate outer-to-inner delegation.
        let inner = Centroided::new(OneSpectrumSource::new(profile_spectrum()));
        let outer = Centroided::new(inner);
        assert_eq!(
            outer.additional_processing_steps(),
            vec![
                CvTerm::from(cv::PEAK_PICKING),
                CvTerm::from(cv::PEAK_PICKING),
            ]
        );
    }

    #[test]
    fn wrapped_source_still_satisfies_conformance_invariants() {
        let mut src = Centroided::new(OneSpectrumSource::new(profile_spectrum()));
        let n = assert_source_invariants(&mut src).expect("conformance");
        assert_eq!(n, 1);
    }

    #[test]
    fn composes_with_the_streaming_mzml_writer() {
        // Centroided is just another SpectrumSource, so it drops straight
        // into write_mzml with no special-casing on the writer's part.
        let mut src = Centroided::new(OneSpectrumSource::new(profile_spectrum()));
        let mut buf = Vec::new();
        crate::write_mzml(&mut src, &mut buf).expect("write_mzml");
        let xml = String::from_utf8(buf).expect("utf8");

        assert!(xml.contains(r#"<spectrumList count="1""#));
        assert!(xml.contains(
            r#"<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>"#
        ));
        assert!(
            !xml.contains(
                r#"<cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>"#
            ),
            "output should not still claim profile mode after centroiding"
        );
        // The dataProcessingList must record that peak picking happened,
        // not just the blanket "Conversion to mzML" step every writer run
        // gets.
        assert!(xml.contains(
            r#"<cvParam cvRef="MS" accession="MS:1000035" name="peak picking" value=""/>"#
        ));
    }

    #[test]
    fn inner_errors_pass_through() {
        let mut inner = OneSpectrumSource::new(profile_spectrum());
        inner.then_fail = true;
        let mut src = Centroided::new(inner);
        let results: Vec<_> = src.iter_spectra().collect();
        assert_eq!(results.len(), 2);
        assert_eq!(
            results[0].as_ref().unwrap().scan_mode,
            Some(ScanMode::Centroid)
        );
        let err = results[1].as_ref().unwrap_err();
        assert!(err.to_string().contains("truncated scan"), "{err}");
    }

    fn centroid_one(rec: SpectrumRecord, leniency: Leniency) -> crate::Result<SpectrumRecord> {
        let mut src = Centroided::new(OneSpectrumSource::new(rec)).with_leniency(leniency);
        let mut out: Vec<_> = src.iter_spectra().collect();
        out.remove(0)
    }

    fn problem(r: crate::Result<SpectrumRecord>) -> RecordProblem {
        match r.unwrap_err() {
            Error::Record { problem, .. } => problem,
            other => panic!("expected Error::Record, got {other:?}"),
        }
    }

    #[test]
    fn unknown_scan_mode_is_an_error_not_assumed_profile() {
        let mut rec = profile_spectrum();
        rec.scan_mode = None;
        assert_eq!(
            problem(centroid_one(rec, Leniency::strict())),
            RecordProblem::UnknownScanMode
        );
    }

    #[test]
    fn unknown_scan_mode_passes_through_marked_when_lenient() {
        let mut rec = profile_spectrum();
        rec.scan_mode = None;
        let original = rec.mz.clone();
        let out = centroid_one(rec, Leniency::strict().with(Check::UnknownScanMode)).unwrap();
        assert_eq!(out.mz, original);
        assert_eq!(out.scan_mode, None);
        assert_eq!(out.warnings.len(), 1);
        assert_eq!(out.warnings[0].kind, WarningKind::CentroidNotApplied);
    }

    #[test]
    fn sparse_peak_list_labelled_profile_is_rejected_not_reweighted() {
        // A centroid peak list mislabelled profile: points are hundreds of
        // m/z apart, so local-maximum weighting would move every peak.
        let mut rec = profile_spectrum();
        rec.mz = vec![150.0, 480.0, 910.0, 1500.0, 2400.0, 3100.0];
        rec.intensity = vec![5.0, 90.0, 20.0, 60.0, 10.0, 40.0];
        let p = problem(centroid_one(rec, Leniency::strict()));
        assert!(matches!(p, RecordProblem::NotProfile { .. }), "{p:?}");
    }

    #[test]
    fn gaps_in_a_zero_stripped_profile_are_not_averaged_across() {
        // Two dense profile bumps with the zeros between them removed, and a
        // lone point far from both. Each result must sit inside its own run.
        let mut rec = profile_spectrum();
        rec.mz = vec![
            500.000, 500.001, 500.002, 500.003, 500.004, // bump A
            600.000, 600.001, 600.002, 600.003, 600.004, // bump B
            700.000, // isolated point
        ];
        rec.intensity = vec![1.0, 5.0, 10.0, 5.0, 1.0, 1.0, 4.0, 8.0, 4.0, 1.0, 3.0];
        let out = centroid_one(rec, Leniency::strict()).unwrap();
        assert_eq!(out.mz.len(), 3, "{:?}", out.mz);
        assert!((out.mz[0] - 500.002).abs() < 0.0015);
        assert!((out.mz[1] - 600.002).abs() < 0.0015);
        assert_eq!(out.mz[2], 700.0);
        assert_eq!(out.intensity, vec![10.0, 8.0, 3.0]);
    }

    #[test]
    fn mismatched_array_lengths_are_an_error_not_truncated() {
        let mut rec = profile_spectrum();
        rec.intensity.push(1.0);
        assert!(matches!(
            problem(centroid_one(rec, Leniency::strict())),
            RecordProblem::ArrayLengthMismatch { .. }
        ));
    }
}
