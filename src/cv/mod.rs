//! Typed PSI-MS and UO controlled-vocabulary terms used by the mzML writer.
//!
//! Every constant here is generated from `psi-ms.obo` by
//! `scripts/gen-cv.py` (the term list is `scripts/cv-terms.txt`), and the
//! generated file records the OBO source URL and `data-version`
//! ([`PSI_MS_VERSION`], [`PSI_MS_SOURCE_URL`]). A [`Term`] can only be built
//! inside this module, so an accession is always paired with its exact OBO
//! name.
//!
//! Reader crates can turn a [`Term`] into a [`CvTerm`](crate::CvTerm) with
//! `From`:
//!
//! ```
//! use openmassspec_core::{cv, CvTerm};
//!
//! let term: CvTerm = cv::PEAK_PICKING.into();
//! assert_eq!(term.accession, "MS:1000035");
//! assert_eq!(term.name, "peak picking");
//! ```

mod terms;

pub use terms::*;

/// A PSI-MS (`MS:`) or Unit Ontology (`UO:`) term known at compile time.
///
/// The fields are private: the only values are the generated constants in
/// this module, so the accession and name always match the ontology.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Term {
    accession: &'static str,
    name: &'static str,
}

impl Term {
    /// Accession, e.g. `"MS:1000016"`.
    pub const fn accession(self) -> &'static str {
        self.accession
    }

    /// Exact OBO name, e.g. `"scan start time"`.
    pub const fn name(self) -> &'static str {
        self.name
    }

    /// mzML `cvRef` for this term: `"UO"` for Unit Ontology accessions,
    /// `"MS"` otherwise.
    pub fn cv_ref(self) -> &'static str {
        if self.accession.starts_with("UO:") {
            "UO"
        } else {
            "MS"
        }
    }
}

/// Look up a generated term by accession.
pub fn lookup(accession: &str) -> Option<Term> {
    ALL.iter().copied().find(|t| t.accession == accession)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lookup_finds_generated_terms() {
        assert_eq!(lookup("MS:1000016"), Some(SCAN_START_TIME));
        assert_eq!(lookup("MS:0000000"), None);
        assert_eq!(MINUTE.cv_ref(), "UO");
        assert_eq!(MZ.cv_ref(), "MS");
    }

    #[test]
    fn accessions_are_unique() {
        let mut seen = std::collections::HashSet::new();
        for t in ALL {
            assert!(seen.insert(t.accession()), "duplicate {}", t.accession());
        }
    }
}
