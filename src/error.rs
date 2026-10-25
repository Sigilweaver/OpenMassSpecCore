//! Aggregate error type for the OpenMassSpec stack.
//!
//! Each vendor crate defines its own narrow `Error` enum. Code that handles
//! errors from several vendors uniformly (the umbrella `openmassspec-io`
//! crate, the `vendor2mzml` CLI, the mzML writer in this crate) works with
//! [`Error`]: a [`SpectrumSource`](crate::SpectrumSource) yields
//! `Result<_, Error>` items, and vendor errors travel inside
//! [`Error::Vendor`].

use std::io;

/// Stack-wide error type.
///
/// Vendor errors are erased into a boxed trait object so this enum does not
/// change as vendor crates evolve. The enum is `#[non_exhaustive]`; match it
/// with a wildcard arm.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// I/O failure (file not found, permission denied, short read, ...).
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    /// Vendor parser failed. The inner error is the vendor crate's own error
    /// type, preserved via `Box<dyn std::error::Error>` so callers can
    /// downcast if they need vendor-specific context.
    #[error("vendor parser error: {0}")]
    Vendor(Box<dyn std::error::Error + Send + Sync + 'static>),

    /// Format detection failed or the input does not match a supported vendor.
    #[error("unsupported or unrecognized format: {0}")]
    Format(String),

    /// Conformance harness rejected a record stream.
    #[error("conformance violation: {0}")]
    Conformance(String),

    /// A spectrum or chromatogram record failed the writer's checks. The
    /// record is identified by its native id (or chromatogram id).
    #[error("invalid record {record_id:?}: {problem}")]
    Record {
        /// `native_id` of the spectrum, or `id` of the chromatogram.
        record_id: String,
        /// What was wrong.
        problem: RecordProblem,
    },
}

/// What is wrong with a record rejected by the mzML writer or the
/// centroider (see [`Error::Record`]).
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum RecordProblem {
    /// Two parallel arrays have different lengths.
    #[error("{array} array has {found} values, expected {expected} (length of {reference})")]
    ArrayLengthMismatch {
        /// The array whose length is wrong.
        array: &'static str,
        /// The array it must match.
        reference: &'static str,
        /// Required length.
        expected: usize,
        /// Length found.
        found: usize,
    },
    /// A float field or array value is NaN or infinite.
    #[error("{field} is not finite ({value})")]
    NonFinite {
        /// Field name, with `[i]` for an array element.
        field: String,
        /// The offending value.
        value: f64,
    },
    /// The m/z array is not in ascending order.
    #[error("m/z array is not ascending at index {index} ({previous} then {value})")]
    UnsortedMz {
        /// Index of the first value that is smaller than its predecessor.
        index: usize,
        /// Value at `index - 1`.
        previous: f64,
        /// Value at `index`.
        value: f64,
    },
    /// The native id was already written for an earlier record.
    #[error("native id already used by an earlier record")]
    DuplicateId,
    /// A precursor `spectrumRef` names a spectrum that is not in the output.
    #[error("precursor references spectrum {reference:?}, which is not in the output")]
    DanglingPrecursorRef {
        /// The referenced native id.
        reference: String,
    },
    /// The reader did not say whether the spectrum is profile or centroid.
    #[error("scan mode is unknown (neither profile nor centroid was declared)")]
    UnknownScanMode,
    /// A spectrum labelled profile is a sparse peak list, not profile data.
    #[error(
        "labelled profile but looks like a peak list \
         (median m/z spacing {median_spacing:.4} is {relative_ppm:.0} ppm of the m/z)"
    )]
    NotProfile {
        /// Median distance between consecutive m/z values.
        median_spacing: f64,
        /// That distance relative to the median m/z, in parts per million.
        relative_ppm: f64,
    },
    /// Charge, intensity or CCS is set but the selected-ion m/z is not.
    #[error("precursor has {field} but no selected-ion m/z")]
    PrecursorWithoutSelectedMz {
        /// Which field is set.
        field: &'static str,
    },
    /// Per-peak mobility needs a mobility kind that the writer cannot express.
    #[error("{0}")]
    Unsupported(String),
}

/// Convenience alias.
pub type Result<T> = std::result::Result<T, Error>;

impl Error {
    /// Wrap a vendor-crate error into [`Error::Vendor`].
    ///
    /// Usable directly as a `map_err` argument:
    /// `decode(scan).map_err(Error::vendor)`.
    pub fn vendor<E>(err: E) -> Self
    where
        E: std::error::Error + Send + Sync + 'static,
    {
        Error::Vendor(Box::new(err))
    }
}
