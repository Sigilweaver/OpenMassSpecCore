//! Vendor-neutral spectrum / chromatogram / run records.
//!
//! These are the canonical handshake between a vendor parser (which knows
//! how to decode its native format) and any downstream consumer (mzML
//! writer, column-store ingest, Python bindings, ...). A vendor parser
//! implements [`crate::SpectrumSource`] which yields these records; a
//! consumer takes any `SpectrumSource` without caring which vendor produced
//! it.

use crate::enums::{Activation, Analyzer, MobilityArrayKind, Polarity, ScanMode};
use std::collections::BTreeMap;

/// A PSI-MS controlled-vocabulary term.
///
/// `accession` is a stable identifier (`MS:NNNNNNN` or `UO:NNNNNNN`); `name`
/// is the human-readable term. mzML output writes these verbatim.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct CvTerm {
    pub accession: &'static str,
    pub name: String,
}

impl<'de> serde::Deserialize<'de> for CvTerm {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        #[derive(serde::Deserialize)]
        struct OwnedCvTerm {
            accession: String,
            name: String,
        }

        let value = OwnedCvTerm::deserialize(deserializer)?;
        Ok(Self {
            accession: intern_accession(value.accession),
            name: value.name,
        })
    }
}

fn intern_accession(value: String) -> &'static str {
    Box::leak(value.into_boxed_str())
}

impl CvTerm {
    pub fn new(accession: &'static str, name: impl Into<String>) -> Self {
        Self {
            accession,
            name: name.into(),
        }
    }
}

/// Precursor metadata for an MS2+ spectrum.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct PrecursorInfo {
    /// Isolation-window center m/z (target).
    pub target_mz: Option<f64>,
    /// Monoisotopic-resolved precursor m/z (selected ion).
    pub selected_mz: Option<f64>,
    /// Total isolation width in m/z. The mzML writer splits this into
    /// symmetric lower/upper offsets around `target_mz`.
    pub isolation_width: Option<f64>,
    pub charge: Option<i32>,
    /// Precursor intensity (when known).
    pub intensity: Option<f64>,
    /// Collision energy. Interpretation depends on `ce_is_nce`.
    pub collision_energy: Option<f64>,
    /// `true` if `collision_energy` is normalized (NCE %); `false` for eV.
    pub ce_is_nce: bool,
    /// Native ID of the precursor scan when known. Vendors should populate
    /// this with the same string they will use as the `native_id` for that
    /// MS1 spectrum, so mzML `spectrumRef` lookups round-trip cleanly.
    pub precursor_native_id: Option<String>,
    pub activation: Option<Activation>,
    /// Analyzer that recorded the precursor scan; needed by mzML to
    /// disambiguate CID vs beam-type CID on FTMS instruments.
    pub analyzer: Option<Analyzer>,
    /// Collision cross-sectional area of the selected ion, in square
    /// angstroms. Derived from ion-mobility measurements (Bruker 1/K0,
    /// Waters drift time, etc); vendors that only expose raw mobility
    /// should convert before populating this rather than leaving it for
    /// downstream consumers to guess a calibration.
    pub ccs: Option<f64>,
}

/// One fully-decoded spectrum.
///
/// Retention time is stored in **seconds** (the mzML preferred unit). Vendors
/// that natively store minutes (e.g. Thermo) convert at the boundary.
///
/// `mz` is f64 and `intensity` is f32 because that is what the PSI-MS CV
/// `64-bit float` / `32-bit float` defaults are and what every downstream
/// search engine expects. Vendors that decode lower-precision arrays should
/// widen to these types.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SpectrumRecord {
    /// Zero-based position in the source file (used as mzML `index=`).
    pub index: usize,
    /// One-based, source-stable scan number. For Bruker bundles this is a
    /// running counter assigned by the iterator (since PASEF frames produce
    /// many spectra per frame).
    pub scan_number: u32,
    /// Verbatim mzML native ID for this spectrum. Vendors populate this with
    /// the appropriate `controllerType=...`, `frame=... scan=...`, or
    /// `function=... process=... scan=...` literal.
    pub native_id: String,
    pub ms_level: u32,
    pub polarity: Option<Polarity>,
    pub scan_mode: Option<ScanMode>,
    pub analyzer: Option<Analyzer>,
    /// Numeric acquisition-event or grouping identifier decoded from the
    /// source file, when available.
    #[serde(default)]
    pub acquisition_event_id: Option<u32>,
    /// Thermo-style scan filter (or vendor-equivalent). Optional; populated
    /// by parsers that have a meaningful filter string.
    pub filter: Option<String>,
    /// Retention time in seconds.
    pub retention_time_sec: f64,
    /// Total ion current. If `None`, the mzML writer computes it from
    /// `intensity`.
    pub total_ion_current: Option<f64>,
    /// Base-peak m/z. If `None`, the mzML writer computes it from `mz` /
    /// `intensity`.
    pub base_peak_mz: Option<f64>,
    /// Base-peak intensity. If `None`, the mzML writer computes it from
    /// `intensity`.
    pub base_peak_intensity: Option<f64>,
    /// Lowest observed m/z. If `None`, the writer uses `mz.first()`.
    pub low_mz: Option<f64>,
    /// Highest observed m/z. If `None`, the writer uses `mz.last()`.
    pub high_mz: Option<f64>,
    pub ion_injection_time_ms: Option<f64>,
    /// Mean inverse reduced ion mobility (1/K0) for the spectrum, when
    /// applicable (Bruker timsTOF, Waters TWIMS).
    pub inv_mobility: Option<f64>,
    /// FAIMS compensation voltage in volts, when the instrument has a FAIMS
    /// Pro interface and reports it per-scan (Thermo Orbitrap Fusion/Lumos/
    /// Eclipse/Ascend). `None` when not applicable or not decoded.
    pub faims_cv: Option<f64>,
    pub precursor: Option<PrecursorInfo>,
    pub mz: Vec<f64>,
    pub intensity: Vec<f32>,
    /// Per-peak inverse reduced ion mobility, parallel to `mz` / `intensity`,
    /// when an IMS-resolved parser opts to preserve it. Length must equal
    /// `mz.len()` when present.
    pub inv_mobility_per_peak: Option<Vec<f32>>,
    /// Namespaced reader-specific values. Keys use `<reader>.<field>`.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub extra: BTreeMap<String, String>,
}

impl SpectrumRecord {
    /// Total ion current, computed from `intensity` when not pre-populated.
    pub fn effective_tic(&self) -> f64 {
        self.total_ion_current
            .unwrap_or_else(|| self.intensity.iter().map(|&v| v as f64).sum())
    }

    /// Base-peak `(mz, intensity)`, computed from the arrays when not
    /// pre-populated. Returns `None` when the spectrum has no peaks.
    pub fn effective_base_peak(&self) -> Option<(f64, f64)> {
        if let (Some(mz), Some(i)) = (self.base_peak_mz, self.base_peak_intensity) {
            return Some((mz, i));
        }
        if self.intensity.is_empty() {
            return None;
        }
        let mut best_idx = 0usize;
        let mut best = self.intensity[0];
        for (i, &v) in self.intensity.iter().enumerate().skip(1) {
            if v > best {
                best = v;
                best_idx = i;
            }
        }
        Some((self.mz[best_idx], best as f64))
    }

    /// `(low_mz, high_mz)`, falling back to the first/last entries in `mz`
    /// when not pre-populated. Returns `None` for empty spectra.
    pub fn effective_mz_range(&self) -> Option<(f64, f64)> {
        match (self.low_mz, self.high_mz) {
            (Some(lo), Some(hi)) => Some((lo, hi)),
            _ => {
                if self.mz.is_empty() {
                    None
                } else {
                    Some((*self.mz.first().unwrap(), *self.mz.last().unwrap()))
                }
            }
        }
    }
}

/// A chromatogram trace (TIC, BPC, SRM/MRM transition).
#[derive(Debug, Clone)]
pub struct ChromatogramRecord {
    pub index: usize,
    pub id: String,
    /// e.g. `"total ion current chromatogram"`, `"basepeak chromatogram"`,
    /// `"selected reaction monitoring chromatogram"`.
    pub chromatogram_type: Option<CvTerm>,
    pub precursor_mz: Option<f64>,
    pub product_mz: Option<f64>,
    /// Retention time in seconds.
    pub time_sec: Vec<f32>,
    pub intensity: Vec<f32>,
}

/// Run-level metadata.
///
/// Vendors construct this once per file. The mzML writer uses it to populate
/// `<fileDescription>`, `<sourceFileList>`, `<instrumentConfigurationList>`
/// and `<softwareList>`.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RunMetadata {
    /// The source file name to put in `<sourceFile name="...">`. Typically
    /// `Path::file_name()`; for directory-based formats (Bruker `.d/`,
    /// Waters `.raw/`), the directory name.
    pub source_file_name: String,
    /// PSI-MS CV term for the source file format, e.g.
    /// `("MS:1000563", "Thermo RAW format")`.
    pub source_file_format: CvTerm,
    /// PSI-MS CV term for the native ID format used in `native_id` fields,
    /// e.g. `("MS:1000768", "Thermo nativeID format")`.
    pub native_id_format: CvTerm,
    /// Instrument CV term. Vendors resolve this from their own model lookup.
    pub instrument: CvTerm,
    /// Serial number of the physical instrument, when the vendor's metadata
    /// exposes one. Emitted as `MS:1000529` ("instrument serial number") on
    /// the default `instrumentConfiguration`.
    pub instrument_serial_number: Option<String>,
    /// Software identifier (e.g. `"opentfraw"`).
    pub software_name: String,
    pub software_version: String,
    /// Name of the original vendor acquisition software that recorded this
    /// run (e.g. `"Xcalibur"`, `"MassHunter"`, `"Analyst"`, `"LabSolutions"`),
    /// when the vendor's own instrument-method metadata exposes one. This is
    /// distinct from `software_name`/`software_version` above, which
    /// identify the vendor-reader crate producing this mzML, not the
    /// instrument-control software that acquired the data. `None` when the
    /// vendor parser does not decode this information. Emitted as a second
    /// `<software>` entry in `<softwareList>` when present.
    pub acquisition_software_name: Option<String>,
    /// Version string for `acquisition_software_name`, when known. Ignored
    /// if `acquisition_software_name` is `None`.
    pub acquisition_software_version: Option<String>,
    /// Instrument acquisition start time (RFC 3339), when available.
    pub start_timestamp: Option<String>,
    /// Interpretation of any per-peak ion mobility array carried by the
    /// spectra in this run. `None` is treated as the Bruker convention
    /// ([`MobilityArrayKind::InverseReducedVsPerCm2`]) for back-compat.
    pub mobility_array_kind: Option<MobilityArrayKind>,
    /// Physical mass-analyzer components this instrument has, in a stable
    /// order. Static instrument configuration, not derived from scan data -
    /// vendors already resolve `instrument` from a model lookup, so this is
    /// the same kind of fixed lookup.
    ///
    /// Each entry gets its own `<instrumentConfiguration>` (in addition to
    /// the default `IC1` built from `instrument`), and any spectrum whose
    /// [`SpectrumRecord::analyzer`](crate::types::SpectrumRecord::analyzer)
    /// matches one of these gets an `instrumentConfigurationRef` pointing at
    /// it, so hybrid instruments (e.g. quad/ion-trap MS2 vs Orbitrap MS1 on
    /// a Fusion) preserve per-scan analyzer identity. Leave empty to keep
    /// the previous single-`IC1`-for-everything behavior.
    pub analyzers: Vec<Analyzer>,
    /// Namespaced reader-specific values. Keys use `<reader>.<field>`.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub extra: BTreeMap<String, String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spectrum() -> SpectrumRecord {
        SpectrumRecord {
            index: 0,
            scan_number: 1,
            native_id: "scan=1".into(),
            ms_level: 1,
            polarity: Some(Polarity::Positive),
            scan_mode: Some(ScanMode::Centroid),
            analyzer: Some(Analyzer::TOFMS),
            acquisition_event_id: Some(17),
            filter: None,
            retention_time_sec: 1.5,
            total_ion_current: None,
            base_peak_mz: None,
            base_peak_intensity: None,
            low_mz: None,
            high_mz: None,
            ion_injection_time_ms: None,
            inv_mobility: None,
            faims_cv: None,
            precursor: None,
            mz: vec![100.0],
            intensity: vec![12.0],
            inv_mobility_per_peak: None,
            extra: BTreeMap::from([("openszraw.event_label".into(), "survey".into())]),
        }
    }

    fn run_metadata() -> RunMetadata {
        RunMetadata {
            source_file_name: "test.raw".into(),
            source_file_format: CvTerm::new("MS:1000563", "Thermo RAW format"),
            native_id_format: CvTerm::new("MS:1000768", "Thermo nativeID format"),
            instrument: CvTerm::new("MS:1001911", "Q Exactive"),
            instrument_serial_number: None,
            software_name: "test-reader".into(),
            software_version: "1.0".into(),
            acquisition_software_name: None,
            acquisition_software_version: None,
            start_timestamp: None,
            mobility_array_kind: None,
            analyzers: Vec::new(),
            extra: BTreeMap::from([("openwraw.source_type".into(), "ESI".into())]),
        }
    }

    #[test]
    fn serde_round_trips_vendor_extras_and_acquisition_event_id() {
        let spectrum: SpectrumRecord =
            serde_json::from_str(&serde_json::to_string(&spectrum()).unwrap()).unwrap();
        assert_eq!(spectrum.acquisition_event_id, Some(17));
        assert_eq!(spectrum.extra["openszraw.event_label"], "survey");

        let metadata: RunMetadata =
            serde_json::from_str(&serde_json::to_string(&run_metadata()).unwrap()).unwrap();
        assert_eq!(metadata.extra["openwraw.source_type"], "ESI");
    }

    #[test]
    fn serde_defaults_new_fields_when_reading_old_records() {
        let mut spectrum_json = serde_json::to_value(spectrum()).unwrap();
        let spectrum_object = spectrum_json.as_object_mut().unwrap();
        spectrum_object.remove("acquisition_event_id");
        spectrum_object.remove("extra");
        let spectrum: SpectrumRecord = serde_json::from_value(spectrum_json).unwrap();
        assert_eq!(spectrum.acquisition_event_id, None);
        assert!(spectrum.extra.is_empty());

        let mut metadata_json = serde_json::to_value(run_metadata()).unwrap();
        metadata_json.as_object_mut().unwrap().remove("extra");
        let metadata: RunMetadata = serde_json::from_value(metadata_json).unwrap();
        assert!(metadata.extra.is_empty());
    }
}
