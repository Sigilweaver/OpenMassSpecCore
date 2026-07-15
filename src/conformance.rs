//! Conformance checks for any [`crate::SpectrumSource`].
//!
//! This module hosts a small invariant suite that every vendor crate runs
//! in its integration tests against a real fixture. It is *not* a
//! validator for arbitrary mzML files; it operates on the in-memory
//! `SpectrumRecord` stream so failures pinpoint the parser, not the
//! writer.
//!
//! Invariants checked (per spectrum unless noted):
//!
//! 1. `mz.len() == intensity.len()`. Empty spectra are allowed because
//!    sparse PASEF / SRM frames can legitimately yield zero peaks within
//!    a scan range.
//! 2. `inv_mobility_per_peak`, when present, has the same length as `mz`.
//! 3. `total_ion_current`, when provided by the parser, equals
//!    `sum(intensity)` within a relative tolerance.
//! 4. `base_peak_intensity`, when provided, equals the maximum of
//!    `intensity` within tolerance.
//! 5. MS2+ spectra carry a [`crate::PrecursorInfo`] with at least one of
//!    `target_mz`, `selected_mz`, or `precursor_native_id` populated.
//! 6. Retention time is non-decreasing within the spectrum stream of one
//!    function / acquisition. The check is per native-ID prefix so
//!    interleaved Waters functions or Bruker MS1/MS2 frames are not
//!    flagged.
//! 7. The first spectrum's `index` is 0 and indices are strictly
//!    increasing.
//!
//! On failure each check emits a [`ConformanceError`] with the offending
//! `native_id` so the operator can jump straight to the problem.

use std::collections::HashMap;

use crate::source::SpectrumSource;
use crate::types::SpectrumRecord;

/// Relative tolerance applied to TIC / base-peak floating-point checks.
const FLOAT_REL_TOL: f64 = 1e-4;

/// Failure modes detected by [`assert_source_invariants`].
#[derive(Debug)]
pub enum ConformanceError {
    /// `mz` / `intensity` arrays have mismatched lengths.
    PeakArrayLengthMismatch {
        native_id: String,
        mz_len: usize,
        intensity_len: usize,
    },
    /// `inv_mobility_per_peak.len()` did not match `mz.len()`.
    MobilityArrayLengthMismatch {
        native_id: String,
        mz_len: usize,
        mobility_len: usize,
    },
    /// `total_ion_current` did not match `sum(intensity)` within tolerance.
    TicMismatch {
        native_id: String,
        declared: f64,
        computed: f64,
    },
    /// `base_peak_intensity` did not match `max(intensity)` within tolerance.
    BasePeakIntensityMismatch {
        native_id: String,
        declared: f64,
        computed: f64,
    },
    /// MS2+ spectrum was missing precursor info.
    MissingPrecursor { native_id: String, ms_level: u32 },
    /// Retention time went backwards within one acquisition stream.
    RetentionTimeNonMonotonic {
        prefix: String,
        previous: f64,
        current: f64,
        native_id: String,
    },
    /// Spectrum index sequence was not strictly increasing or did not start at 0.
    IndexSequence {
        native_id: String,
        previous: Option<usize>,
        current: usize,
    },
    /// Spectrum had no peaks at all.
    EmptySpectrum { native_id: String },
}

impl std::fmt::Display for ConformanceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::PeakArrayLengthMismatch {
                native_id,
                mz_len,
                intensity_len,
            } => write!(
                f,
                "{native_id}: mz.len()={mz_len} != intensity.len()={intensity_len}"
            ),
            Self::MobilityArrayLengthMismatch {
                native_id,
                mz_len,
                mobility_len,
            } => write!(
                f,
                "{native_id}: inv_mobility_per_peak.len()={mobility_len} != mz.len()={mz_len}"
            ),
            Self::TicMismatch {
                native_id,
                declared,
                computed,
            } => {
                write!(
                    f,
                    "{native_id}: declared TIC {declared} != computed {computed}"
                )
            }
            Self::BasePeakIntensityMismatch {
                native_id,
                declared,
                computed,
            } => write!(
                f,
                "{native_id}: declared base-peak intensity {declared} != max(intensity) {computed}"
            ),
            Self::MissingPrecursor {
                native_id,
                ms_level,
            } => write!(f, "{native_id}: ms_level={ms_level} but no precursor info"),
            Self::RetentionTimeNonMonotonic {
                prefix,
                previous,
                current,
                native_id,
            } => write!(
                f,
                "{native_id}: RT regressed in stream '{prefix}' ({previous} -> {current})"
            ),
            Self::IndexSequence {
                native_id,
                previous,
                current,
            } => write!(
                f,
                "{native_id}: index sequence broken (prev={previous:?}, current={current})"
            ),
            Self::EmptySpectrum { native_id } => write!(f, "{native_id}: empty spectrum"),
        }
    }
}

impl std::error::Error for ConformanceError {}

fn rel_close(a: f64, b: f64, tol: f64) -> bool {
    let scale = a.abs().max(b.abs()).max(1.0);
    (a - b).abs() <= tol * scale
}

/// Per-native-ID-prefix key used to scope retention-time monotonicity.
///
/// We split on the first whitespace and strip any trailing `scan=...`
/// token so that mzML native IDs like `function=2 process=0 scan=17`,
/// `frame=12 scan=8`, and `controllerType=0 controllerNumber=1 scan=42`
/// each group correctly.
fn rt_stream_key(native_id: &str) -> String {
    native_id
        .split_whitespace()
        .filter(|tok| !tok.starts_with("scan="))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Run the invariant suite against a borrowed source. Returns the first
/// failure encountered, or `Ok(spectrum_count)` if every spectrum passed.
///
/// The source is iterated to completion on success; on failure iteration
/// stops at the offending spectrum.
pub fn assert_source_invariants<S: SpectrumSource>(src: &mut S) -> Result<usize, ConformanceError> {
    assert_iter_invariants(src.iter_spectra())
}

/// Run the invariant suite against any iterator yielding spectra. Useful
/// for vendors whose top-level mzML writer does not yet implement
/// [`SpectrumSource`].
pub fn assert_iter_invariants<I: IntoIterator<Item = SpectrumRecord>>(
    iter: I,
) -> Result<usize, ConformanceError> {
    let mut last_index: Option<usize> = None;
    let mut last_rt: HashMap<String, f64> = HashMap::new();
    let mut count = 0usize;
    for spectrum in iter {
        check_one(&spectrum, last_index, &mut last_rt)?;
        last_index = Some(spectrum.index);
        count += 1;
    }
    Ok(count)
}

fn check_one(
    s: &SpectrumRecord,
    last_index: Option<usize>,
    last_rt: &mut HashMap<String, f64>,
) -> Result<(), ConformanceError> {
    if s.mz.len() != s.intensity.len() {
        return Err(ConformanceError::PeakArrayLengthMismatch {
            native_id: s.native_id.clone(),
            mz_len: s.mz.len(),
            intensity_len: s.intensity.len(),
        });
    }
    if let Some(mob) = &s.inv_mobility_per_peak {
        if mob.len() != s.mz.len() {
            return Err(ConformanceError::MobilityArrayLengthMismatch {
                native_id: s.native_id.clone(),
                mz_len: s.mz.len(),
                mobility_len: mob.len(),
            });
        }
    }
    if let Some(tic) = s.total_ion_current {
        let computed: f64 = s.intensity.iter().map(|&v| v as f64).sum();
        if !rel_close(tic, computed, FLOAT_REL_TOL) {
            return Err(ConformanceError::TicMismatch {
                native_id: s.native_id.clone(),
                declared: tic,
                computed,
            });
        }
    }
    if let Some(bp) = s.base_peak_intensity {
        if s.intensity.is_empty() {
            // Nothing to compare against; trust the parser.
            let _ = bp;
        } else {
            let computed = s
                .intensity
                .iter()
                .copied()
                .fold(f32::NEG_INFINITY, f32::max) as f64;
            if !rel_close(bp, computed, FLOAT_REL_TOL) {
                return Err(ConformanceError::BasePeakIntensityMismatch {
                    native_id: s.native_id.clone(),
                    declared: bp,
                    computed,
                });
            }
        }
    }
    if s.ms_level >= 2 {
        let has_precursor = match &s.precursor {
            Some(p) => {
                p.target_mz.is_some() || p.selected_mz.is_some() || p.precursor_native_id.is_some()
            }
            None => false,
        };
        if !has_precursor {
            return Err(ConformanceError::MissingPrecursor {
                native_id: s.native_id.clone(),
                ms_level: s.ms_level,
            });
        }
    }
    match last_index {
        None => {
            if s.index != 0 {
                return Err(ConformanceError::IndexSequence {
                    native_id: s.native_id.clone(),
                    previous: None,
                    current: s.index,
                });
            }
        }
        Some(prev) => {
            if s.index <= prev {
                return Err(ConformanceError::IndexSequence {
                    native_id: s.native_id.clone(),
                    previous: Some(prev),
                    current: s.index,
                });
            }
        }
    }
    let key = rt_stream_key(&s.native_id);
    if let Some(prev) = last_rt.get(&key).copied() {
        // Allow tiny floating-point regressions (<1 us).
        if s.retention_time_sec + 1e-6 < prev {
            return Err(ConformanceError::RetentionTimeNonMonotonic {
                prefix: key,
                previous: prev,
                current: s.retention_time_sec,
                native_id: s.native_id.clone(),
            });
        }
    }
    last_rt.insert(key, s.retention_time_sec);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::RunMetadata;
    use crate::{Polarity, PrecursorInfo, SpectrumRecord};

    struct ToySource(Vec<SpectrumRecord>);

    impl SpectrumSource for ToySource {
        fn run_metadata(&self) -> RunMetadata {
            RunMetadata {
                source_file_name: "toy".into(),
                source_file_format: crate::CvTerm {
                    accession: "MS:1000563",
                    name: "Thermo RAW format".into(),
                },
                native_id_format: crate::CvTerm {
                    accession: "MS:1000768",
                    name: "Thermo nativeID format".into(),
                },
                instrument: crate::CvTerm {
                    accession: "MS:1000483",
                    name: "Thermo Fisher Scientific instrument model".into(),
                },
                software_name: "test".into(),
                software_version: "0.0".into(),
                start_timestamp: None,
                mobility_array_kind: None,
            }
        }
        fn iter_spectra<'a>(&'a mut self) -> Box<dyn Iterator<Item = SpectrumRecord> + 'a> {
            Box::new(self.0.clone().into_iter())
        }
    }

    fn ok_spec(index: usize, ms_level: u32, rt: f64) -> SpectrumRecord {
        let mz = vec![100.0, 200.0];
        let intensity = vec![10.0f32, 20.0];
        SpectrumRecord {
            index,
            scan_number: (index + 1) as u32,
            native_id: format!("controllerType=0 controllerNumber=1 scan={}", index + 1),
            ms_level,
            polarity: Some(Polarity::Positive),
            scan_mode: None,
            analyzer: None,
            filter: None,
            retention_time_sec: rt,
            total_ion_current: Some(30.0),
            base_peak_mz: Some(200.0),
            base_peak_intensity: Some(20.0),
            low_mz: Some(100.0),
            high_mz: Some(200.0),
            ion_injection_time_ms: None,
            inv_mobility: None,
            faims_cv: None,
            precursor: if ms_level >= 2 {
                Some(PrecursorInfo {
                    selected_mz: Some(500.0),
                    ..Default::default()
                })
            } else {
                None
            },
            mz,
            intensity,
            inv_mobility_per_peak: None,
        }
    }

    #[test]
    fn happy_path_passes() {
        let mut src = ToySource(vec![
            ok_spec(0, 1, 1.0),
            ok_spec(1, 2, 2.0),
            ok_spec(2, 1, 3.0),
        ]);
        assert_eq!(assert_source_invariants(&mut src).unwrap(), 3);
    }

    #[test]
    fn detects_tic_mismatch() {
        let mut s = ok_spec(0, 1, 1.0);
        s.total_ion_current = Some(999.0);
        let mut src = ToySource(vec![s]);
        let err = assert_source_invariants(&mut src).unwrap_err();
        assert!(matches!(err, ConformanceError::TicMismatch { .. }));
    }

    #[test]
    fn detects_missing_precursor() {
        let mut s = ok_spec(0, 2, 1.0);
        s.precursor = None;
        let mut src = ToySource(vec![s]);
        let err = assert_source_invariants(&mut src).unwrap_err();
        assert!(matches!(err, ConformanceError::MissingPrecursor { .. }));
    }

    #[test]
    fn detects_rt_regression_within_same_stream() {
        let s0 = ok_spec(0, 1, 5.0);
        let s1 = ok_spec(1, 1, 2.0);
        let mut src = ToySource(vec![s0, s1]);
        let err = assert_source_invariants(&mut src).unwrap_err();
        assert!(matches!(
            err,
            ConformanceError::RetentionTimeNonMonotonic { .. }
        ));
    }

    #[test]
    fn detects_mobility_length_mismatch() {
        let mut s = ok_spec(0, 1, 1.0);
        s.inv_mobility_per_peak = Some(vec![0.5]);
        let mut src = ToySource(vec![s]);
        let err = assert_source_invariants(&mut src).unwrap_err();
        assert!(matches!(
            err,
            ConformanceError::MobilityArrayLengthMismatch { .. }
        ));
    }
}
