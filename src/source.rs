//! The `SpectrumSource` trait: every vendor parser implements this.

use crate::types::{ChromatogramRecord, RunMetadata, SpectrumRecord};

/// A source of decoded mass spectra.
///
/// Vendors implement this on whatever value carries their open file state
/// (e.g. `RawFileReader` + a `&mut Read+Seek` source for opentfraw, a
/// `Reader` for opentimstdf).
///
/// The trait deliberately uses boxed iterators rather than RPITIT so that
/// implementations can pick a different underlying iterator type per call
/// without leaking that into the trait signature, and so consumers can hold
/// a `&mut dyn SpectrumSource` for downstream plumbing (mzML writer, ingest
/// pipelines, language bindings).
pub trait SpectrumSource {
    /// Run-level metadata. Cheap to call; vendors typically build this once.
    fn run_metadata(&self) -> RunMetadata;

    /// Iterate every spectrum the file contains. Spectra the parser cannot
    /// decode should be skipped silently; the writer trusts whatever the
    /// iterator yields.
    ///
    /// The iterator borrows `self` mutably so vendors can stream from disk
    /// without buffering the whole run in memory.
    fn iter_spectra<'a>(&'a mut self) -> Box<dyn Iterator<Item = SpectrumRecord> + 'a>;

    /// Iterate chromatogram traces (TIC, BPC, SRM). Defaults to an empty
    /// iterator; most parsers do not synthesize chromatograms.
    fn iter_chromatograms<'a>(&'a mut self) -> Box<dyn Iterator<Item = ChromatogramRecord> + 'a> {
        Box::new(std::iter::empty())
    }

    /// Total number of spectra the source will yield, when known cheaply.
    /// Used by the mzML writer to populate `<spectrumList count="...">`. If
    /// `None`, the writer falls back to buffering spectrum offsets and
    /// patching the count at the end.
    fn spectrum_count_hint(&self) -> Option<usize> {
        None
    }
}
