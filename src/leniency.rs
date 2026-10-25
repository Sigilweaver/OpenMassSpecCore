//! Opt-in lenient handling of inputs that strict mode rejects.
//!
//! By default the writer and the centroider return an error when a record is
//! not what they expect. A [`Leniency`] names the specific [`Check`]s that
//! may instead continue, marking each affected record with a
//! [`Warning`](crate::Warning). An empty set is strict.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::fmt;
use std::str::FromStr;

/// One strictness check that lenient mode can relax.
///
/// The built-in checks are variants. A reader crate defines its own with
/// [`Check::other`] and a namespaced code (`"opensxraw.calibration_range"`),
/// so a new check needs no Core release. Built-in names (`as_str`,
/// `FromStr`) are what command lines and Python keywords accept; a reader
/// parses its own names itself.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum Check {
    /// A spectrum whose scan mode is unknown: written with no
    /// profile/centroid term, and left unchanged by centroiding. The
    /// resulting mzML lacks the spectrum-representation term mzML expects,
    /// so it can fail semantic validation by design.
    UnknownScanMode,
    /// A precursor that names a spectrum not in the output: the
    /// `spectrumRef` is dropped. References are resolved against spectra
    /// written earlier, so a reference to a later spectrum counts as
    /// dangling.
    DanglingPrecursorRef,
    /// A precursor with a collision energy and no unit: the value is
    /// written as a unit-less userParam instead of a cvParam.
    UnknownCollisionEnergyUnit,
    /// The source announced more spectra than it yielded and did not report
    /// an error: the write succeeds with a warning on the run instead of
    /// failing.
    SpectrumCountShortfall,
    /// A check defined by a reader crate, identified by a namespaced code.
    Other(Cow<'static, str>),
}

impl Check {
    /// Every built-in check this version knows. A newer Core can add
    /// checks, which makes [`Leniency::all`] more lenient.
    pub const ALL: &'static [Check] = &[
        Check::UnknownScanMode,
        Check::DanglingPrecursorRef,
        Check::UnknownCollisionEnergyUnit,
        Check::SpectrumCountShortfall,
    ];

    /// A reader-defined check with a namespaced `code`. A code that names a
    /// built-in check gives that built-in variant, so `allows` never misses
    /// it.
    pub fn other(code: impl Into<Cow<'static, str>>) -> Self {
        let code = code.into();
        match code.parse::<Check>() {
            Ok(builtin) => builtin,
            Err(_) => Check::Other(code),
        }
    }

    /// Stable name (`"unknown_scan_mode"`, ...); the code of an
    /// [`Other`](Check::Other) check.
    pub fn as_str(&self) -> &str {
        match self {
            Self::UnknownScanMode => "unknown_scan_mode",
            Self::DanglingPrecursorRef => "dangling_precursor_ref",
            Self::UnknownCollisionEnergyUnit => "unknown_collision_energy_unit",
            Self::SpectrumCountShortfall => "spectrum_count_shortfall",
            Self::Other(code) => code,
        }
    }
}

impl fmt::Display for Check {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// A check name that matches no built-in [`Check`].
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("unknown check {name:?}; known checks: {known}")]
#[non_exhaustive]
pub struct UnknownCheck {
    /// The name that was given.
    pub name: String,
    /// Comma-separated list of the valid names.
    pub known: String,
}

impl FromStr for Check {
    type Err = UnknownCheck;

    /// Parse a built-in check name. Reader-defined names are rejected here
    /// on purpose, so a typo fails; readers build [`Check::other`] from the
    /// names they own.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Check::ALL
            .iter()
            .find(|c| c.as_str() == s)
            .cloned()
            .ok_or_else(|| UnknownCheck {
                name: s.to_string(),
                known: Check::ALL
                    .iter()
                    .map(|c| c.as_str())
                    .collect::<Vec<_>>()
                    .join(", "),
            })
    }
}

/// The set of [`Check`]s that may be relaxed. Empty (the default) is strict.
///
/// ```
/// use openmassspec_core::{Check, Leniency};
///
/// let lenient = Leniency::strict().with(Check::UnknownScanMode);
/// assert!(lenient.allows(&Check::UnknownScanMode));
/// assert!(!lenient.allows(&Check::DanglingPrecursorRef));
/// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct Leniency {
    allowed: BTreeSet<Check>,
}

impl Leniency {
    /// No check is relaxed.
    pub fn strict() -> Self {
        Self::default()
    }

    /// Every built-in check this version knows is relaxed. Reader-defined
    /// ([`Check::Other`]) checks are not included: add them with
    /// [`with`](Self::with).
    pub fn all() -> Self {
        Self {
            allowed: Check::ALL.iter().cloned().collect(),
        }
    }

    /// Relax `check`.
    #[must_use]
    pub fn with(mut self, check: Check) -> Self {
        self.allowed.insert(check);
        self
    }

    /// `true` if `check` may be relaxed.
    pub fn allows(&self, check: &Check) -> bool {
        self.allowed.contains(check)
    }

    /// `true` if no check is relaxed.
    pub fn is_strict(&self) -> bool {
        self.allowed.is_empty()
    }

    /// The relaxed checks, in a stable order.
    pub fn iter(&self) -> impl Iterator<Item = &Check> + '_ {
        self.allowed.iter()
    }
}

impl FromIterator<Check> for Leniency {
    fn from_iter<I: IntoIterator<Item = Check>>(iter: I) -> Self {
        Self {
            allowed: iter.into_iter().collect(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_round_trip_and_unknown_names_list_the_valid_ones() {
        for c in Check::ALL {
            assert_eq!(&c.as_str().parse::<Check>().unwrap(), c);
        }
        let e = "nope".parse::<Check>().unwrap_err();
        assert!(e.to_string().contains("unknown_scan_mode"));
    }

    #[test]
    fn strict_by_default() {
        assert!(Leniency::default().is_strict());
        assert!(!Leniency::all().is_strict());
        assert!(Leniency::all().allows(&Check::DanglingPrecursorRef));
        let mine = Check::other("opensxraw.calibration_range");
        assert!(!Leniency::all().allows(&mine));
        assert!(Leniency::strict().with(mine.clone()).allows(&mine));
        assert!("opensxraw.calibration_range".parse::<Check>().is_err());
        // A built-in code gives the built-in variant.
        assert_eq!(Check::other("unknown_scan_mode"), Check::UnknownScanMode);
        assert!(Leniency::strict()
            .with(Check::UnknownScanMode)
            .allows(&Check::other("unknown_scan_mode")));
    }
}
