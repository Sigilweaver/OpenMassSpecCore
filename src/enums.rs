//! Shared enums and types used across the open mass-spec parsers.
//!
//! These mirror the enums historically defined in `opentfraw::types` (which is
//! where the vocabulary first lived) but are intentionally **vendor neutral**:
//! anything Thermo-specific (e.g. firmware byte codes, scan-filter symbols)
//! stays in the vendor crate.

/// Detector / mass analyzer family.
///
/// The variants follow Thermo's preamble encoding (the most detailed source
/// available across vendors); other vendors map their detector class onto the
/// closest equivalent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Analyzer {
    /// Ion trap.
    ITMS,
    /// Triple quadrupole.
    TQMS,
    /// Single quadrupole.
    SQMS,
    /// Time-of-flight.
    TOFMS,
    /// Fourier-transform (Orbitrap, FT-ICR).
    FTMS,
    /// Magnetic sector.
    Sector,
}

/// Scan polarity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Polarity {
    Negative,
    Positive,
}

/// Whether a spectrum holds centroided peak picks or the raw profile signal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum ScanMode {
    Centroid,
    Profile,
}

/// MSn order (MS1, MS2, ...).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MsPower {
    Undefined,
    Ms1,
    Ms2,
    Ms3,
    Ms4,
    Ms5,
    Ms6,
    Ms7,
    Ms8,
}

/// Activation method used for MSn fragmentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Activation {
    /// Higher-energy collisional dissociation (Q Exactive style).
    HCD,
    /// Multi-photon induced dissociation.
    MPID,
    /// Electron transfer dissociation.
    ETD,
    /// Collision-induced dissociation. On FTMS analyzers this is
    /// conventionally rendered as beam-type CID (HCD-equivalent) in mzML.
    CID,
    /// Electron-capture dissociation.
    ECD,
    /// Infrared multiphoton dissociation.
    IRMPD,
    /// Proton-transfer / activated-ion variant.
    PD,
    /// Pulsed q dissociation.
    PQD,
    /// Ultraviolet photodissociation.
    UVPD,
    /// Surface-induced dissociation.
    SID,
    /// ETD with supplemental HCD.
    EThcD,
}

/// Compression codec applied to a `<binaryDataArray>`'s `<binary>` payload
/// by the mzML writer.
///
/// Selects both the bytes actually written (compressed or not) and the
/// matching PSI-MS cvParam. `Zlib` is the default for
/// [`crate::write_mzml`]/[`crate::write_indexed_mzml`], matching
/// msconvert/OpenMS output; pass `NoCompression` to
/// [`crate::write_mzml_with_compression`]/
/// [`crate::write_indexed_mzml_with_compression`] to keep the previous
/// (pre-1.4) uncompressed byte layout.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Compression {
    /// No compression: `MS:1000576` ("no compression").
    NoCompression,
    /// DEFLATE/zlib compression: `MS:1000574` ("zlib compression").
    #[default]
    Zlib,
}

/// Unit + meaning of a per-peak ion mobility array.
///
/// Selects the CV term and unit emitted by the mzML writer when a
/// [`crate::SpectrumRecord::inv_mobility_per_peak`] array is present.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum MobilityArrayKind {
    /// Bruker TIMS convention: per-peak 1/K0 in volt-second per square
    /// centimeter. Emitted as MS:1003008 "raw inverse reduced ion mobility
    /// array" with unit MS:1002814.
    InverseReducedVsPerCm2,
    /// Waters traveling-wave IMS convention: per-peak drift time in
    /// milliseconds. Emitted as MS:1003007 "raw ion mobility array" with
    /// unit UO:0000028 "millisecond".
    DriftTimeMilliseconds,
}

impl MsPower {
    /// Numeric MS level (1 for MS1, 2 for MS2, ...). Returns 1 for `Undefined`,
    /// matching the convention used in mzML output.
    pub fn ms_level(self) -> u32 {
        match self {
            Self::Undefined | Self::Ms1 => 1,
            Self::Ms2 => 2,
            Self::Ms3 => 3,
            Self::Ms4 => 4,
            Self::Ms5 => 5,
            Self::Ms6 => 6,
            Self::Ms7 => 7,
            Self::Ms8 => 8,
        }
    }
}
