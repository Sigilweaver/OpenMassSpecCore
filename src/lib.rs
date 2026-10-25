//! `openmassspec-core` is the shared foundation for the open Rust mass-spec
//! parsers (`opentfraw`, `opentimstdf`, `openwraw`, `openaraw`, `opensxraw`,
//! `openszraw`).
//!
//! It exposes:
//!
//! * Vendor-neutral spectrum / chromatogram / run records:
//!   [`SpectrumRecord`], [`PrecursorInfo`], [`ChromatogramRecord`],
//!   [`RunMetadata`], [`CvTerm`]. All are `#[non_exhaustive]` and built with
//!   constructors plus `with_*` setters.
//! * Shared enums: [`Polarity`], [`ScanMode`], [`Analyzer`], [`MsPower`],
//!   [`Activation`], [`Compression`], [`MobilityArrayKind`].
//! * A trait every vendor parser implements: [`SpectrumSource`]. Its
//!   iterators yield `Result` items, so decode errors reach the consumer.
//! * One canonical mzML 1.1.0 writer: [`write_mzml`] and
//!   [`write_indexed_mzml`] (zlib-compressed binary data arrays, stop at the
//!   first record error), with [`WriteOptions`] for codec and
//!   [`RecordErrorPolicy`] control. Each call returns a [`WriteSummary`].
//! * Typed PSI-MS / UO terms generated from `psi-ms.obo`: [`cv`].
//! * A streaming profile-to-centroid transform, [`Centroided`], that wraps
//!   any `SpectrumSource`.
//! * A conformance harness for parser test suites: [`conformance`].
//!
//! The Apache Arrow bridge lives in the separate `openmassspec-arrow` crate.
//!
//! Each vendor crate is a standalone tool (a user can pull in `opentfraw`
//! alone and get parsing **and** mzML export); `openmassspec-core` is the
//! shared vocabulary that keeps the parsers in lock-step.

#![warn(missing_docs)]

mod centroid;
mod enums;
mod error;
mod leniency;
mod mzml;
mod source;
mod types;
mod validate;
mod warning;

pub mod conformance;
pub mod cv;

pub use centroid::Centroided;
pub use enums::{
    Activation, Analyzer, CollisionEnergyUnit, Compression, MobilityArrayKind, MsPower,
    ParseLabelError, Polarity, ScanMode,
};
pub use error::{Error, RecordProblem, Result};
pub use leniency::{Check, Leniency, UnknownCheck};
pub use mzml::{
    write_indexed_mzml, write_indexed_mzml_with_options, write_mzml, write_mzml_with_options,
    RecordErrorPolicy, WriteIssue, WriteOptions, WriteSummary,
};
pub use source::{ChromatogramIter, SpectrumIter, SpectrumSource};
pub use types::{ChromatogramRecord, CvTerm, PrecursorInfo, RunMetadata, SpectrumRecord};
pub use warning::{Warning, WarningKind};

/// Crate version (set from `Cargo.toml`).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
