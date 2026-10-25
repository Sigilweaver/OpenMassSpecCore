//! Shared enums used across the open mass-spec parsers.
//!
//! The enums are **vendor neutral**: anything vendor-specific (e.g. firmware
//! byte codes, scan-filter symbols) stays in the vendor crate.
//!
//! Every enum is `#[non_exhaustive]`, so new variants can be added in a minor
//! release. Code outside this crate that matches on one needs a wildcard arm.

/// Detector / mass analyzer family.
///
/// The variants follow Thermo's preamble encoding (the most detailed source
/// available across vendors); other vendors map their detector class onto the
/// closest equivalent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
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
#[non_exhaustive]
pub enum Polarity {
    /// Negative ion mode.
    Negative,
    /// Positive ion mode.
    Positive,
}

/// Whether a spectrum holds centroided peak picks or the raw profile signal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum ScanMode {
    /// Centroided peak list.
    Centroid,
    /// Continuum (profile) signal.
    Profile,
}

/// MSn order (MS1, MS2, ...).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum MsPower {
    /// MS order not recorded.
    Undefined,
    /// MS1.
    Ms1,
    /// MS2.
    Ms2,
    /// MS3.
    Ms3,
    /// MS4.
    Ms4,
    /// MS5.
    Ms5,
    /// MS6.
    Ms6,
    /// MS7.
    Ms7,
    /// MS8.
    Ms8,
}

/// Activation method used for MSn fragmentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum Activation {
    /// Higher-energy collisional dissociation (Q Exactive style).
    HCD,
    /// Multi-photon induced dissociation.
    MPID,
    /// Electron transfer dissociation.
    ETD,
    /// Collision-induced dissociation (resonant). Beam-type CID is
    /// [`Activation::HCD`]; the analyzer is not used to choose.
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
/// matching PSI-MS cvParam. `Zlib` is the default; select another codec with
/// [`WriteOptions::with_compression`](crate::WriteOptions::with_compression).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
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
#[non_exhaustive]
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

impl Polarity {
    /// Stable lowercase label (`"positive"`, `"negative"`), used as the
    /// `polarity` value in `openmassspec-arrow` batches.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Positive => "positive",
            Self::Negative => "negative",
        }
    }
}

impl ScanMode {
    /// Stable lowercase label (`"centroid"`, `"profile"`), used as the
    /// `scan_mode` value in `openmassspec-arrow` batches.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Centroid => "centroid",
            Self::Profile => "profile",
        }
    }
}

impl Analyzer {
    /// Stable lowercase label (`"itms"`, `"tqms"`, `"sqms"`, `"tof"`,
    /// `"ftms"`, `"sector"`), used as the `analyzer` value in
    /// `openmassspec-arrow` batches.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ITMS => "itms",
            Self::TQMS => "tqms",
            Self::SQMS => "sqms",
            Self::TOFMS => "tof",
            Self::FTMS => "ftms",
            Self::Sector => "sector",
        }
    }
}

impl Activation {
    /// Stable lowercase label (`"hcd"`, `"cid"`, `"ethcd"`, ...), used as the
    /// `precursor_activation` value in `openmassspec-arrow` batches.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::CID => "cid",
            Self::HCD => "hcd",
            Self::ETD => "etd",
            Self::ECD => "ecd",
            Self::UVPD => "uvpd",
            Self::PQD => "pqd",
            Self::PD => "pd",
            Self::SID => "sid",
            Self::EThcD => "ethcd",
            Self::IRMPD => "irmpd",
            Self::MPID => "mpid",
        }
    }
}

impl MobilityArrayKind {
    /// Stable label (`"inverse_reduced_k0"`, `"drift_time_ms"`), used as the
    /// `mobility_array_kind` value in `openmassspec-arrow` batches.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::InverseReducedVsPerCm2 => "inverse_reduced_k0",
            Self::DriftTimeMilliseconds => "drift_time_ms",
        }
    }
}
