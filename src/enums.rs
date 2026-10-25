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
    #[serde(rename = "itms")]
    ITMS,
    /// Triple quadrupole.
    #[serde(rename = "tqms")]
    TQMS,
    /// Single quadrupole.
    #[serde(rename = "sqms")]
    SQMS,
    /// Time-of-flight.
    #[serde(rename = "tof")]
    TOFMS,
    /// Fourier-transform (Orbitrap, FT-ICR). The mzML writer cannot tell
    /// them apart and writes the Orbitrap term, so an FT-ICR run is labelled
    /// Orbitrap until a variant that names the analyzer exists.
    #[serde(rename = "ftms")]
    FTMS,
    /// Magnetic sector.
    #[serde(rename = "sector")]
    Sector,
}

/// Scan polarity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum Polarity {
    /// Negative ion mode.
    #[serde(rename = "negative")]
    Negative,
    /// Positive ion mode.
    #[serde(rename = "positive")]
    Positive,
}

/// Whether a spectrum holds centroided peak picks or the raw profile signal.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum ScanMode {
    /// Centroided peak list.
    #[serde(rename = "centroid")]
    Centroid,
    /// Continuum (profile) signal.
    #[serde(rename = "profile")]
    Profile,
}

/// Activation method used for MSn fragmentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum Activation {
    /// Higher-energy collisional dissociation (Q Exactive style).
    #[serde(rename = "hcd")]
    HCD,
    /// Multi-photon induced dissociation.
    #[serde(rename = "mpid")]
    MPID,
    /// Electron transfer dissociation.
    #[serde(rename = "etd")]
    ETD,
    /// Collision-induced dissociation (resonant). Beam-type CID is
    /// [`Activation::HCD`]; the analyzer is not used to choose.
    #[serde(rename = "cid")]
    CID,
    /// Electron-capture dissociation.
    #[serde(rename = "ecd")]
    ECD,
    /// Infrared multiphoton dissociation.
    #[serde(rename = "irmpd")]
    IRMPD,
    /// Proton-transfer / activated-ion variant.
    #[serde(rename = "pd")]
    PD,
    /// Pulsed q dissociation.
    #[serde(rename = "pqd")]
    PQD,
    /// Ultraviolet photodissociation.
    #[serde(rename = "uvpd")]
    UVPD,
    /// Surface-induced dissociation.
    #[serde(rename = "sid")]
    SID,
    /// ETD with supplemental HCD.
    #[serde(rename = "ethcd")]
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
    #[serde(rename = "inverse_reduced_k0")]
    InverseReducedVsPerCm2,
    /// Waters traveling-wave IMS convention: per-peak drift time in
    /// milliseconds. Emitted as MS:1003007 "raw ion mobility array" with
    /// unit UO:0000028 "millisecond".
    #[serde(rename = "drift_time_ms")]
    DriftTimeMilliseconds,
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

/// Unit of [`PrecursorInfo::collision_energy`](crate::PrecursorInfo::collision_energy).
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub enum CollisionEnergyUnit {
    /// Electronvolts (absolute collision energy).
    #[serde(rename = "ev")]
    Ev,
    /// Normalized collision energy, percent.
    #[serde(rename = "nce")]
    Nce,
}

impl CollisionEnergyUnit {
    /// Stable label (`"ev"`, `"nce"`).
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Ev => "ev",
            Self::Nce => "nce",
        }
    }
}

/// A label that is not the `as_str` label of any variant of the named enum.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("unknown {enum_name} label {label:?}")]
#[non_exhaustive]
pub struct ParseLabelError {
    /// Name of the enum that was being parsed.
    pub enum_name: &'static str,
    /// The label that matched no variant.
    pub label: String,
}

macro_rules! from_str_labels {
    ($ty:ident: $($variant:ident),+ $(,)?) => {
        impl std::str::FromStr for $ty {
            type Err = ParseLabelError;

            /// Parse the label written by `as_str` (and by serde).
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                $(if s == Self::$variant.as_str() {
                    return Ok(Self::$variant);
                })+
                Err(ParseLabelError {
                    enum_name: stringify!($ty),
                    label: s.to_string(),
                })
            }
        }
    };
}

from_str_labels!(Polarity: Positive, Negative);
from_str_labels!(ScanMode: Centroid, Profile);
from_str_labels!(Analyzer: ITMS, TQMS, SQMS, TOFMS, FTMS, Sector);
from_str_labels!(
    Activation: CID, HCD, ETD, ECD, UVPD, PQD, PD, SID, EThcD, IRMPD, MPID
);
from_str_labels!(MobilityArrayKind: InverseReducedVsPerCm2, DriftTimeMilliseconds);
from_str_labels!(CollisionEnergyUnit: Ev, Nce);

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip<T>(variants: &[T])
    where
        T: Copy
            + PartialEq
            + std::fmt::Debug
            + std::str::FromStr
            + serde::Serialize
            + serde::de::DeserializeOwned,
        <T as std::str::FromStr>::Err: std::fmt::Debug,
    {
        for &v in variants {
            let json = serde_json::to_string(&v).unwrap();
            let label: String = serde_json::from_str(&json).unwrap();
            assert_eq!(label.parse::<T>().unwrap(), v, "{json}");
            assert_eq!(serde_json::from_str::<T>(&json).unwrap(), v);
        }
    }

    #[test]
    fn serde_labels_are_the_as_str_labels() {
        use Activation::*;
        assert_eq!(serde_json::to_string(&Analyzer::TOFMS).unwrap(), r#""tof""#);
        assert_eq!(
            serde_json::to_string(&Polarity::Positive).unwrap(),
            r#""positive""#
        );
        assert_eq!(
            serde_json::to_string(&MobilityArrayKind::InverseReducedVsPerCm2).unwrap(),
            r#""inverse_reduced_k0""#
        );
        roundtrip(&[Polarity::Positive, Polarity::Negative]);
        roundtrip(&[ScanMode::Centroid, ScanMode::Profile]);
        roundtrip(&[
            Analyzer::ITMS,
            Analyzer::TQMS,
            Analyzer::SQMS,
            Analyzer::TOFMS,
            Analyzer::FTMS,
            Analyzer::Sector,
        ]);
        roundtrip(&[CID, HCD, ETD, ECD, UVPD, PQD, PD, SID, EThcD, IRMPD, MPID]);
        roundtrip(&[
            MobilityArrayKind::InverseReducedVsPerCm2,
            MobilityArrayKind::DriftTimeMilliseconds,
        ]);
        roundtrip(&[CollisionEnergyUnit::Ev, CollisionEnergyUnit::Nce]);
    }

    #[test]
    fn unknown_labels_are_errors() {
        let e = "TOFMS".parse::<Analyzer>().unwrap_err();
        assert_eq!(e.enum_name, "Analyzer");
        assert!(serde_json::from_str::<Analyzer>(r#""TOFMS""#).is_err());
    }
}
