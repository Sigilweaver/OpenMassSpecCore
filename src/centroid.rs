//! Vendor-neutral profile-to-centroid peak picking.
//!
//! Centroiding is a transform over the arrays a [`SpectrumSource`] already
//! yields (`mz`/`intensity`, and optionally `inv_mobility_per_peak`) - once a
//! [`SpectrumRecord`] exists, the operation is identical regardless of which
//! vendor parser produced it. [`Centroided`] wraps any `SpectrumSource` and
//! applies the transform lazily, one spectrum at a time, so it composes with
//! the streaming mzML writer and the Arrow bridge without buffering a whole
//! run in memory.

use crate::enums::ScanMode;
use crate::source::SpectrumSource;
use crate::types::{ChromatogramRecord, RunMetadata, SpectrumRecord};

/// Wraps a [`SpectrumSource`], centroiding every profile-mode spectrum it
/// yields. Spectra already tagged [`ScanMode::Centroid`] pass through
/// unchanged (idempotent).
pub struct Centroided<S: SpectrumSource> {
    inner: S,
    min_intensity: f32,
}

impl<S: SpectrumSource> Centroided<S> {
    /// Wrap `inner`, picking every local-maximum peak regardless of height.
    pub fn new(inner: S) -> Self {
        Self {
            inner,
            min_intensity: 0.0,
        }
    }

    /// Discard picked peaks below `min_intensity` (a simple noise floor).
    pub fn with_min_intensity(mut self, min_intensity: f32) -> Self {
        self.min_intensity = min_intensity;
        self
    }
}

impl<S: SpectrumSource> SpectrumSource for Centroided<S> {
    fn run_metadata(&self) -> RunMetadata {
        self.inner.run_metadata()
    }

    fn iter_spectra<'a>(&'a mut self) -> Box<dyn Iterator<Item = SpectrumRecord> + 'a> {
        let min_intensity = self.min_intensity;
        Box::new(
            self.inner
                .iter_spectra()
                .map(move |rec| centroid_record(rec, min_intensity)),
        )
    }

    fn iter_chromatograms<'a>(&'a mut self) -> Box<dyn Iterator<Item = ChromatogramRecord> + 'a> {
        self.inner.iter_chromatograms()
    }

    fn spectrum_count_hint(&self) -> Option<usize> {
        self.inner.spectrum_count_hint()
    }
}

fn centroid_record(mut rec: SpectrumRecord, min_intensity: f32) -> SpectrumRecord {
    if rec.scan_mode == Some(ScanMode::Centroid) {
        return rec;
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
    rec
}

/// Local-maxima peak picking: a point is a picked peak if it is no smaller
/// than both neighbors and strictly larger than at least one of them (this
/// also correctly picks a single-point spectrum, `n == 1`, since the "beats
/// a neighbor" check is vacuously satisfied). The picked m/z (and, when
/// present, inverse mobility) is the intensity-weighted centroid over the
/// apex and its immediate neighbors; the picked intensity is the apex
/// height. This is intentionally simple - see `timsrust-centroid` /
/// `pyteomics` for more sophisticated peer approaches if this ever needs to
/// improve on plain local-maxima picking.
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

    for i in 0..n {
        let no_smaller_than_left = i == 0 || intensity[i] >= intensity[i - 1];
        let no_smaller_than_right = i == n - 1 || intensity[i] >= intensity[i + 1];
        let beats_a_neighbor = n == 1
            || (i > 0 && intensity[i] > intensity[i - 1])
            || (i < n - 1 && intensity[i] > intensity[i + 1]);
        if !(no_smaller_than_left && no_smaller_than_right && beats_a_neighbor) {
            continue;
        }
        if intensity[i] < min_intensity {
            continue;
        }

        let lo = i.saturating_sub(1);
        let hi = (i + 1).min(n - 1);
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
    use crate::{CvTerm, MsPower, Polarity};

    struct OneSpectrumSource {
        meta: RunMetadata,
        rec: Option<SpectrumRecord>,
    }

    impl OneSpectrumSource {
        fn new(rec: SpectrumRecord) -> Self {
            Self {
                meta: RunMetadata {
                    source_file_name: "toy.raw".into(),
                    source_file_format: CvTerm::new("MS:1000563", "Thermo RAW format"),
                    native_id_format: CvTerm::new("MS:1000768", "Thermo nativeID format"),
                    instrument: CvTerm::new("MS:1001911", "Q Exactive"),
                    software_name: "toy-writer".into(),
                    software_version: "0.0.0".into(),
                    start_timestamp: None,
                    mobility_array_kind: None,
                },
                rec: Some(rec),
            }
        }
    }

    impl SpectrumSource for OneSpectrumSource {
        fn run_metadata(&self) -> RunMetadata {
            self.meta.clone()
        }

        fn iter_spectra<'a>(&'a mut self) -> Box<dyn Iterator<Item = SpectrumRecord> + 'a> {
            Box::new(self.rec.take().into_iter())
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
            filter: None,
            retention_time_sec: 1.0,
            total_ion_current: Some(intensity.iter().map(|&v| v as f64).sum()),
            base_peak_mz: Some(100.0),
            base_peak_intensity: Some(100.0),
            low_mz: Some(99.0),
            high_mz: Some(201.0),
            ion_injection_time_ms: None,
            inv_mobility: None,
            faims_cv: None,
            precursor: None,
            mz,
            intensity,
            inv_mobility_per_peak: None,
        }
    }

    #[test]
    fn centroids_profile_spectrum_at_known_apexes() {
        let mut src = Centroided::new(OneSpectrumSource::new(profile_spectrum()));
        let recs: Vec<_> = src.iter_spectra().collect();
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
        let recs: Vec<_> = src.iter_spectra().collect();
        assert_eq!(recs[0].mz, original_mz);
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
    }
}
