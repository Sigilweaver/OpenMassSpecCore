//! `openmassspec-core` is the shared foundation for the open Rust mass-spec
//! parsers (`opentfraw`, `opentimstdf`, `openwraw`).
//!
//! It exposes:
//!
//! * Vendor-neutral spectrum / chromatogram / run records:
//!   [`SpectrumRecord`], [`PrecursorInfo`], [`ChromatogramRecord`],
//!   [`RunMetadata`], [`CvTerm`].
//! * Shared enums: [`Polarity`], [`ScanMode`], [`Analyzer`], [`MsPower`],
//!   [`Activation`], [`Compression`].
//! * A trait every vendor parser implements: [`SpectrumSource`].
//! * One canonical mzML 1.1.0 writer: [`write_mzml`] and
//!   [`write_indexed_mzml`] (zlib-compressed binary data arrays by default;
//!   [`write_mzml_with_compression`]/[`write_indexed_mzml_with_compression`]
//!   for explicit codec control).
//! * A streaming profile-to-centroid transform, [`Centroided`], that wraps
//!   any `SpectrumSource`.
//!
//! Each vendor crate is a standalone tool (a user can pull in `opentfraw`
//! alone and get parsing **and** mzML export); `openmassspec-core` is the
//! shared vocabulary that keeps the three parsers in lock-step.

mod centroid;
mod enums;
mod error;
mod mzml;
mod source;
mod types;

#[cfg(feature = "arrow")]
pub mod arrow;
pub mod conformance;

pub use centroid::Centroided;
pub use enums::{
    Activation, Analyzer, Compression, MobilityArrayKind, MsPower, Polarity, ScanMode,
};
pub use error::{Error, Result};
pub use mzml::{
    write_indexed_mzml, write_indexed_mzml_with_compression, write_mzml,
    write_mzml_with_compression,
};
pub use source::SpectrumSource;
pub use types::{ChromatogramRecord, CvTerm, PrecursorInfo, RunMetadata, SpectrumRecord};

/// Crate version (set from `Cargo.toml`).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
