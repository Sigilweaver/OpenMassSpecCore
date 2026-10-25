//! The `SpectrumSource` trait: every vendor parser implements this.

use crate::error::Result;
use crate::types::{ChromatogramRecord, CvTerm, RunMetadata, SpectrumRecord};

/// Boxed, fallible spectrum iterator returned by
/// [`SpectrumSource::iter_spectra`]. It is `Send`, so a source can be
/// driven from a worker thread (the Python wheels do).
pub type SpectrumIter<'a> = Box<dyn Iterator<Item = Result<SpectrumRecord>> + Send + 'a>;

/// Boxed, fallible chromatogram iterator returned by
/// [`SpectrumSource::iter_chromatograms`]. It is `Send`.
pub type ChromatogramIter<'a> = Box<dyn Iterator<Item = Result<ChromatogramRecord>> + Send + 'a>;

/// A source of decoded mass spectra.
///
/// Vendors implement this on whatever value carries their open file state
/// (e.g. `RawFileReader` + a `&mut Read+Seek` source for opentfraw, a
/// `Reader` for opentimstdf).
///
/// The trait uses boxed iterators rather than RPITIT so that implementations
/// can pick a different underlying iterator type per call without leaking
/// that into the trait signature, and so consumers can hold a
/// `&mut dyn SpectrumSource` for downstream plumbing (mzML writer, ingest
/// pipelines, language bindings).
///
/// # Errors
///
/// Iterators yield `Result` items with the stack-wide
/// [`Error`](crate::Error). A record the parser cannot decode is yielded as
/// `Err`, never skipped: wrap the vendor error with
/// [`Error::vendor`](crate::Error::vendor) (or a `From` impl in the vendor
/// crate). One error type for every source keeps `dyn SpectrumSource`
/// usable across vendors and lets adapters such as
/// [`Centroided`](crate::Centroided) pass errors through unchanged.
/// Consumers decide what an error means: the mzML writer stops at the first
/// one unless the caller opts into
/// [`RecordErrorPolicy::Skip`](crate::RecordErrorPolicy::Skip).
///
/// # Example
///
/// ```
/// use openmassspec_core::{
///     cv, CvTerm, Error, RunMetadata, SpectrumIter, SpectrumRecord, SpectrumSource,
/// };
///
/// /// Stand-in for a vendor crate's own error type.
/// #[derive(Debug)]
/// struct DecodeError(u32);
/// impl std::fmt::Display for DecodeError {
///     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
///         write!(f, "scan {} is truncated", self.0)
///     }
/// }
/// impl std::error::Error for DecodeError {}
///
/// struct MyReader {
///     scans: Vec<u32>,
/// }
///
/// impl MyReader {
///     fn decode(&self, i: usize, scan: u32) -> Result<SpectrumRecord, DecodeError> {
///         if scan == 0 {
///             return Err(DecodeError(scan));
///         }
///         Ok(SpectrumRecord::new(i, scan, format!("scan={scan}"))
///             .with_ms_level(1)
///             .with_retention_time_sec(f64::from(scan))
///             .with_peaks(vec![100.0, 200.0], vec![1.0, 2.0]))
///     }
/// }
///
/// impl SpectrumSource for MyReader {
///     fn run_metadata(&self) -> RunMetadata {
///         RunMetadata::new(
///             "run.raw",
///             CvTerm::new("MS:1000563", "Thermo RAW format"),
///             CvTerm::new("MS:1000768", "Thermo nativeID format"),
///             CvTerm::new("MS:1000483", "Thermo Fisher Scientific instrument model"),
///             "my-reader",
///             "1.0.0",
///         )
///     }
///
///     fn iter_spectra<'a>(&'a mut self) -> SpectrumIter<'a> {
///         let scans = self.scans.clone();
///         Box::new(
///             scans
///                 .into_iter()
///                 .enumerate()
///                 .map(move |(i, scan)| self.decode(i, scan).map_err(Error::vendor)),
///         )
///     }
/// }
///
/// let mut reader = MyReader { scans: vec![1, 0] };
/// let results: Vec<_> = reader.iter_spectra().collect();
/// assert!(results[0].is_ok());
/// assert!(results[1].is_err());
/// ```
pub trait SpectrumSource {
    /// Run-level metadata. Cheap to call; vendors typically build this once.
    fn run_metadata(&self) -> RunMetadata;

    /// Iterate every spectrum the file contains, in file order.
    ///
    /// A spectrum the parser cannot decode is yielded as `Err` (see the
    /// trait-level docs). Iteration may continue after an `Err` when the
    /// parser can resume at the next spectrum; a parser that cannot resume
    /// ends the iterator after yielding the error.
    ///
    /// The iterator borrows `self` mutably so vendors can stream from disk
    /// without buffering the whole run in memory.
    fn iter_spectra<'a>(&'a mut self) -> SpectrumIter<'a>;

    /// Iterate chromatogram traces (TIC, BPC, SRM), with the same error
    /// contract as [`iter_spectra`](Self::iter_spectra). Defaults to an
    /// empty iterator.
    fn iter_chromatograms<'a>(&'a mut self) -> ChromatogramIter<'a> {
        Box::new(std::iter::empty())
    }

    /// Total number of spectra the source will yield, when known cheaply.
    /// A sizing hint for consumers. The mzML writer does not use it: its
    /// `<spectrumList count>` is the number of spectra actually written.
    fn spectrum_count_hint(&self) -> Option<usize> {
        None
    }

    /// Extra PSI-MS `<processingMethod>` steps, beyond the writer's own
    /// default "Conversion to mzML" step, that belong in the emitted
    /// `<dataProcessingList>`, in the order they should be recorded.
    ///
    /// Defaults to none. Adapters that transform spectra in a way a
    /// consumer needs recorded in the file's processing history (e.g.
    /// [`Centroided`](crate::Centroided) doing peak picking) override this,
    /// delegating to the wrapped source first so steps accumulate through a
    /// chain of adapters. Use the typed constants in [`crate::cv`] where one
    /// exists (`cv::PEAK_PICKING.into()`).
    fn additional_processing_steps(&self) -> Vec<CvTerm> {
        Vec::new()
    }
}
