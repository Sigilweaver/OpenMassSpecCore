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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Polarity {
    Negative,
    Positive,
}

/// Whether a spectrum holds centroided peak picks or the raw profile signal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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
