//! Per-record checks run by the mzML writer and the centroider before a
//! record is used, so a malformed record is an error instead of output that
//! disagrees with itself.

use crate::error::RecordProblem;
use crate::types::{ChromatogramRecord, PrecursorInfo, SpectrumRecord};

fn finite(field: &str, value: f64) -> Result<(), RecordProblem> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(RecordProblem::NonFinite {
            field: field.to_string(),
            value,
        })
    }
}

/// Like [`finite`], but builds the field name only when the check fails.
fn finite_with(field: impl FnOnce() -> String, value: f64) -> Result<(), RecordProblem> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(RecordProblem::NonFinite {
            field: field(),
            value,
        })
    }
}

fn finite_opt(field: &str, value: Option<f64>) -> Result<(), RecordProblem> {
    value.map_or(Ok(()), |v| finite(field, v))
}

/// Reject mismatched array lengths, non-finite numbers and an unsorted m/z
/// array. Does not look at scan mode, ids or analyzers; the writer checks
/// those against its own state.
pub(crate) fn check_spectrum(rec: &SpectrumRecord) -> Result<(), RecordProblem> {
    check_array_lengths(rec)?;
    finite_opt("retention_time_sec", rec.retention_time_sec)?;
    finite_opt("total_ion_current", rec.total_ion_current)?;
    finite_opt("base_peak_mz", rec.base_peak_mz)?;
    finite_opt("base_peak_intensity", rec.base_peak_intensity)?;
    finite_opt("low_mz", rec.low_mz)?;
    finite_opt("high_mz", rec.high_mz)?;
    finite_opt("ion_injection_time_ms", rec.ion_injection_time_ms)?;
    finite_opt("inv_mobility", rec.inv_mobility)?;
    finite_opt("faims_cv", rec.faims_cv)?;
    for (i, w) in rec.scan_windows.iter().enumerate() {
        finite_with(|| format!("scan_windows[{i}].lower"), w.lower)?;
        finite_with(|| format!("scan_windows[{i}].upper"), w.upper)?;
        if w.lower < 0.0 || w.lower > w.upper {
            return Err(RecordProblem::InvalidScanWindow {
                index: i,
                lower: w.lower,
                upper: w.upper,
            });
        }
    }
    if let Some(p) = &rec.precursor {
        check_precursor(p)?;
    }
    for (i, &v) in rec.mz.iter().enumerate() {
        finite(&format!("mz[{i}]"), v)?;
    }
    for (i, &v) in rec.intensity.iter().enumerate() {
        finite(&format!("intensity[{i}]"), f64::from(v))?;
    }
    if let Some(mob) = &rec.inv_mobility_per_peak {
        for (i, &v) in mob.iter().enumerate() {
            finite(&format!("inv_mobility_per_peak[{i}]"), f64::from(v))?;
        }
    }
    for (i, w) in rec.mz.windows(2).enumerate() {
        if w[1] < w[0] {
            return Err(RecordProblem::UnsortedMz {
                index: i + 1,
                previous: w[0],
                value: w[1],
            });
        }
    }
    Ok(())
}

/// Only the array-length checks (used by the centroider, which must not
/// index past the shorter array).
pub(crate) fn check_array_lengths(rec: &SpectrumRecord) -> Result<(), RecordProblem> {
    if rec.intensity.len() != rec.mz.len() {
        return Err(RecordProblem::ArrayLengthMismatch {
            array: "intensity",
            reference: "mz",
            expected: rec.mz.len(),
            found: rec.intensity.len(),
        });
    }
    if let Some(mob) = &rec.inv_mobility_per_peak {
        if mob.len() != rec.mz.len() {
            return Err(RecordProblem::ArrayLengthMismatch {
                array: "inv_mobility_per_peak",
                reference: "mz",
                expected: rec.mz.len(),
                found: mob.len(),
            });
        }
    }
    Ok(())
}

fn check_precursor(p: &PrecursorInfo) -> Result<(), RecordProblem> {
    finite_opt("precursor.target_mz", p.target_mz)?;
    finite_opt("precursor.selected_mz", p.selected_mz)?;
    finite_opt("precursor.isolation_width", p.isolation_width)?;
    finite_opt("precursor.intensity", p.intensity)?;
    for (i, ce) in p.collision_energies.iter().enumerate() {
        finite_with(|| format!("precursor.collision_energies[{i}]"), ce.value)?;
    }
    finite_opt("precursor.ccs", p.ccs)?;
    if p.selected_mz.is_none() {
        for (field, set) in [
            ("charge", p.charge.is_some()),
            ("intensity", p.intensity.is_some()),
            ("ccs", p.ccs.is_some()),
        ] {
            if set {
                return Err(RecordProblem::PrecursorWithoutSelectedMz { field });
            }
        }
    }
    Ok(())
}

pub(crate) fn check_chromatogram(rec: &ChromatogramRecord) -> Result<(), RecordProblem> {
    if rec.intensity.len() != rec.time_sec.len() {
        return Err(RecordProblem::ArrayLengthMismatch {
            array: "intensity",
            reference: "time_sec",
            expected: rec.time_sec.len(),
            found: rec.intensity.len(),
        });
    }
    finite_opt("precursor_mz", rec.precursor_mz)?;
    finite_opt("product_mz", rec.product_mz)?;
    for (i, &v) in rec.time_sec.iter().enumerate() {
        finite(&format!("time_sec[{i}]"), f64::from(v))?;
    }
    for (i, &v) in rec.intensity.iter().enumerate() {
        finite(&format!("intensity[{i}]"), f64::from(v))?;
    }
    Ok(())
}

/// Median distance between consecutive m/z values and that distance as a
/// fraction of the median m/z. `None` for fewer than three points.
pub(crate) fn median_spacing(mz: &[f64]) -> Option<(f64, f64)> {
    if mz.len() < 3 {
        return None;
    }
    let mut gaps: Vec<f64> = mz.windows(2).map(|w| w[1] - w[0]).collect();
    gaps.sort_by(f64::total_cmp);
    let spacing = gaps[gaps.len() / 2];
    let mut sorted = mz.to_vec();
    sorted.sort_by(f64::total_cmp);
    let mid = sorted[sorted.len() / 2];
    if mid <= 0.0 {
        return None;
    }
    Some((spacing, spacing / mid))
}
