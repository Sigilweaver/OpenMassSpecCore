//! Per-record warnings.
//!
//! A [`Warning`] marks a record that was read or written in a lenient mode
//! (see [`Leniency`](crate::Leniency)) or that carries a value the reader is
//! not fully sure of. Warnings travel on the record
//! ([`SpectrumRecord::warnings`](crate::SpectrumRecord::warnings),
//! [`ChromatogramRecord::warnings`](crate::ChromatogramRecord::warnings),
//! [`RunMetadata::warnings`](crate::RunMetadata::warnings)), are written to
//! mzML as `openmassspec: warning <code>` userParams and to Arrow as a
//! `warnings` list column of codes, so a consumer sees them without a log.

use std::borrow::Cow;
use std::fmt;

/// What kind of warning a [`Warning`] is.
///
/// Reader-specific kinds use [`WarningKind::Other`] with a namespaced code
/// such as `"opentimstdf.tof_near_limit"`, so a new kind needs no Core
/// release. Consumers must treat unknown codes as [`WarningKind::Other`].
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum WarningKind {
    /// The reader did not say whether the spectrum is profile or centroid.
    UnknownScanMode,
    /// A precursor `spectrumRef` named a spectrum that is not in the output
    /// and was dropped.
    DanglingPrecursorRef,
    /// Centroiding was requested but not applied to this spectrum.
    CentroidNotApplied,
    /// The precursor has a collision energy whose unit the reader did not
    /// state; it was written without a unit.
    UnknownCollisionEnergyUnit,
    /// The source announced more spectra than were written or skipped.
    SpectrumCountShortfall,
    /// A code defined by a reader crate or a newer Core.
    Other(Cow<'static, str>),
}

impl WarningKind {
    /// Stable machine-readable code (`"unknown_scan_mode"`, ...).
    pub fn code(&self) -> &str {
        match self {
            Self::UnknownScanMode => "unknown_scan_mode",
            Self::DanglingPrecursorRef => "dangling_precursor_ref",
            Self::CentroidNotApplied => "centroid_not_applied",
            Self::UnknownCollisionEnergyUnit => "unknown_collision_energy_unit",
            Self::SpectrumCountShortfall => "spectrum_count_shortfall",
            Self::Other(code) => code,
        }
    }

    /// The kind for `code`; codes this version does not know become
    /// [`WarningKind::Other`], so reading a newer producer's output never
    /// fails.
    pub fn from_code(code: &str) -> Self {
        match code {
            "unknown_scan_mode" => Self::UnknownScanMode,
            "dangling_precursor_ref" => Self::DanglingPrecursorRef,
            "centroid_not_applied" => Self::CentroidNotApplied,
            "unknown_collision_energy_unit" => Self::UnknownCollisionEnergyUnit,
            "spectrum_count_shortfall" => Self::SpectrumCountShortfall,
            other => Self::Other(Cow::Owned(other.to_string())),
        }
    }
}

impl fmt::Display for WarningKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.code())
    }
}

impl serde::Serialize for WarningKind {
    fn serialize<S: serde::Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        s.serialize_str(self.code())
    }
}

impl<'de> serde::Deserialize<'de> for WarningKind {
    fn deserialize<D: serde::Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        let code = String::deserialize(d)?;
        Ok(Self::from_code(&code))
    }
}

/// A warning attached to a record.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct Warning {
    /// What kind of warning this is.
    pub kind: WarningKind,
    /// Human-readable detail.
    pub message: String,
}

impl Warning {
    /// A warning of `kind` with a detail message.
    pub fn new(kind: WarningKind, message: impl Into<String>) -> Self {
        Self {
            kind,
            message: message.into(),
        }
    }

    /// A reader-specific warning with a namespaced `code`.
    pub fn other(code: impl Into<Cow<'static, str>>, message: impl Into<String>) -> Self {
        Self::new(WarningKind::Other(code.into()), message)
    }
}

impl fmt::Display for Warning {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.kind, self.message)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unknown_codes_round_trip_as_other() {
        let k = WarningKind::from_code("opentimstdf.tof_near_limit");
        assert_eq!(
            k,
            WarningKind::Other(Cow::Borrowed("opentimstdf.tof_near_limit"))
        );
        assert_eq!(k.code(), "opentimstdf.tof_near_limit");
        assert_eq!(
            WarningKind::from_code("unknown_scan_mode"),
            WarningKind::UnknownScanMode
        );
    }

    #[test]
    fn serde_uses_codes_and_accepts_unknown_ones() {
        let w = Warning::new(WarningKind::UnknownScanMode, "no mode");
        let json = serde_json::to_string(&w).unwrap();
        assert_eq!(json, r#"{"kind":"unknown_scan_mode","message":"no mode"}"#);
        let back: Warning = serde_json::from_str(&json).unwrap();
        assert_eq!(back, w);
        let future: Warning =
            serde_json::from_str(r#"{"kind":"from_the_future","message":"x"}"#).unwrap();
        assert_eq!(future.kind.code(), "from_the_future");
    }
}
