//! The record checks and the error types are usable from outside the crate.

use openmassspec_core::{ChromatogramRecord, Error, RecordProblem, SpectrumRecord};
use std::error::Error as _;

#[test]
fn validate_reports_typed_problems() {
    let ok = SpectrumRecord::new(0, 1, "scan=1").with_peaks(vec![1.0, 2.0], vec![3.0, 4.0]);
    assert_eq!(ok.validate(), Ok(()));

    let short = SpectrumRecord::new(0, 1, "scan=1").with_peaks(vec![1.0, 2.0], vec![3.0]);
    assert!(matches!(
        short.validate(),
        Err(RecordProblem::ArrayLengthMismatch { .. })
    ));

    let nan = SpectrumRecord::new(0, 1, "scan=1").with_peaks(vec![1.0], vec![f32::NAN]);
    assert!(matches!(
        nan.validate(),
        Err(RecordProblem::NonFinite { .. })
    ));

    let unsorted = SpectrumRecord::new(0, 1, "scan=1").with_peaks(vec![2.0, 1.0], vec![1.0, 1.0]);
    assert!(matches!(
        unsorted.validate(),
        Err(RecordProblem::UnsortedMz { index: 1, .. })
    ));

    let level0 = SpectrumRecord::new(0, 1, "scan=1").with_ms_level(0);
    assert_eq!(level0.validate(), Err(RecordProblem::MsLevelZero));

    let chrom = ChromatogramRecord::new(0, "TIC", vec![0.0, 1.0], vec![1.0]);
    assert!(matches!(
        chrom.validate(),
        Err(RecordProblem::ArrayLengthMismatch { .. })
    ));
}

#[test]
fn readers_can_build_record_errors_without_a_core_release() {
    let problem =
        RecordProblem::other("opentimstdf.mobility_range", "1/K0 outside the calibration");
    let err = Error::record("frame=1 scan=3", problem.clone());
    assert!(err.to_string().contains("opentimstdf.mobility_range"));
    assert!(matches!(err, Error::Record { problem: p, .. } if p == problem));
}

#[test]
fn vendor_errors_keep_their_source_chain() {
    let err = Error::vendor(std::io::Error::other("block 7 is corrupt"));
    let source = err.source().expect("vendor error has a source");
    assert_eq!(source.to_string(), "block 7 is corrupt");
}
