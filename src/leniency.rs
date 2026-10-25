//! Opt-in lenient handling of inputs that strict mode rejects.
//!
//! By default the writer and the centroider return an error when a record is
//! not what they expect. A [`Leniency`] names the specific [`Check`]s that
//! may instead continue, marking each affected record with a
//! [`Warning`](crate::Warning). An empty set is strict.

use std::collections::BTreeSet;
use std::fmt;
use std::str::FromStr;

/// One strictness check that lenient mode can relax.
///
/// New checks are added as variants; their [`as_str`](Check::as_str) names
/// are what command lines and Python keywords accept.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum Check {
    /// A spectrum whose scan mode is unknown: written with no
    /// profile/centroid term, and left unchanged by centroiding.
    UnknownScanMode,
    /// A precursor that names a spectrum not in the output: the
    /// `spectrumRef` is dropped.
    DanglingPrecursorRef,
}

impl Check {
    /// Every check this version knows.
    pub const ALL: &'static [Check] = &[Check::UnknownScanMode, Check::DanglingPrecursorRef];

    /// Stable name (`"unknown_scan_mode"`, ...).
    pub fn as_str(self) -> &'static str {
        match self {
            Self::UnknownScanMode => "unknown_scan_mode",
            Self::DanglingPrecursorRef => "dangling_precursor_ref",
        }
    }
}

impl fmt::Display for Check {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// A check name that matches no [`Check`].
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

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Check::ALL
            .iter()
            .copied()
            .find(|c| c.as_str() == s)
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
/// assert!(lenient.allows(Check::UnknownScanMode));
/// assert!(!lenient.allows(Check::DanglingPrecursorRef));
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

    /// Every check this version knows is relaxed.
    pub fn all() -> Self {
        Self {
            allowed: Check::ALL.iter().copied().collect(),
        }
    }

    /// Relax `check`.
    #[must_use]
    pub fn with(mut self, check: Check) -> Self {
        self.allowed.insert(check);
        self
    }

    /// `true` if `check` may be relaxed.
    pub fn allows(&self, check: Check) -> bool {
        self.allowed.contains(&check)
    }

    /// `true` if no check is relaxed.
    pub fn is_strict(&self) -> bool {
        self.allowed.is_empty()
    }

    /// The relaxed checks, in a stable order.
    pub fn iter(&self) -> impl Iterator<Item = Check> + '_ {
        self.allowed.iter().copied()
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
        for &c in Check::ALL {
            assert_eq!(c.as_str().parse::<Check>().unwrap(), c);
        }
        let e = "nope".parse::<Check>().unwrap_err();
        assert!(e.to_string().contains("unknown_scan_mode"));
    }

    #[test]
    fn strict_by_default() {
        assert!(Leniency::default().is_strict());
        assert!(!Leniency::all().is_strict());
        assert!(Leniency::all().allows(Check::DanglingPrecursorRef));
    }
}
