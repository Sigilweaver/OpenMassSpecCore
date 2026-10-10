//! Vendor-neutral spectrum / chromatogram / run records.
//!
//! These are the canonical handshake between a vendor parser (which knows
//! how to decode its native format) and any downstream consumer (mzML
//! writer, column-store ingest, Python bindings, ...). A vendor parser
//! implements [`crate::SpectrumSource`] which yields these records; a
//! consumer takes any `SpectrumSource` without caring which vendor produced
//! it.
//!
//! Every record struct is `#[non_exhaustive]`, so fields can be added in a
//! minor release. Code outside this crate builds records with the
//! constructors and `with_*` setters (or assigns the public fields of a
//! constructed value); struct literals and `..Default::default()` only work
//! inside this crate.

use std::borrow::Cow;
use std::collections::BTreeMap;

use crate::cv;
use crate::enums::{Activation, Analyzer, MobilityArrayKind, Polarity, ScanMode};

/// Generates builder setters: `KIND SETTER FIELD: TYPE`.
///
/// * `val` takes the field type as is.
/// * `into` takes `impl Into<TYPE>`.
/// * `opt` sets an `Option<TYPE>` field and takes `impl Into<Option<TYPE>>`,
///   so both `T` and `Option<T>` are accepted.
macro_rules! setters {
    ($ty:ty { $( $kind:ident $setter:ident $field:ident : $t:ty ),* $(,)? }) => {
        impl $ty {
            $( setters!(@one $kind $setter $field : $t); )*
        }
    };
    (@one val $setter:ident $field:ident : $t:ty) => {
        #[doc = concat!("Set `", stringify!($field), "`.")]
        #[must_use]
        pub fn $setter(mut self, value: $t) -> Self {
            self.$field = value;
            self
        }
    };
    (@one into $setter:ident $field:ident : $t:ty) => {
        #[doc = concat!("Set `", stringify!($field), "`.")]
        #[must_use]
        pub fn $setter(mut self, value: impl Into<$t>) -> Self {
            self.$field = value.into();
            self
        }
    };
    (@one opt $setter:ident $field:ident : $t:ty) => {
        #[doc = concat!(
            "Set `", stringify!($field), "`. Accepts a value or an `Option`."
        )]
        #[must_use]
        pub fn $setter(mut self, value: impl Into<Option<$t>>) -> Self {
            self.$field = value.into();
            self
        }
    };
}

/// A PSI-MS controlled-vocabulary term.
///
/// `accession` is a stable identifier (`MS:NNNNNNN` or `UO:NNNNNNN`); `name`
/// is the human-readable term. mzML output writes these verbatim. Both are
/// `Cow<'static, str>`: terms built from literals or from the generated
/// [`crate::cv`] constants borrow, deserialized terms own their strings.
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct CvTerm {
    /// Accession, e.g. `"MS:1000563"`.
    pub accession: Cow<'static, str>,
    /// Exact CV name, e.g. `"Thermo RAW format"`.
    pub name: Cow<'static, str>,
}

impl CvTerm {
    /// Build a term from an accession and its exact CV name.
    pub fn new(
        accession: impl Into<Cow<'static, str>>,
        name: impl Into<Cow<'static, str>>,
    ) -> Self {
        Self {
            accession: accession.into(),
            name: name.into(),
        }
    }
}

impl From<cv::Term> for CvTerm {
    fn from(term: cv::Term) -> Self {
        Self::new(term.accession(), term.name())
    }
}

/// Precursor metadata for an MS2+ spectrum.
///
/// Build with `PrecursorInfo::default()` and the `with_*` setters.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct PrecursorInfo {
    /// Isolation-window center m/z (target).
    pub target_mz: Option<f64>,
    /// Monoisotopic-resolved precursor m/z (selected ion).
    pub selected_mz: Option<f64>,
    /// Total isolation width in m/z. The mzML writer splits this into
    /// symmetric lower/upper offsets around `target_mz`.
    pub isolation_width: Option<f64>,
    /// Precursor charge state.
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
    /// Activation (dissociation) method.
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
///
/// # Example
///
/// A reader in another crate builds a record with [`SpectrumRecord::new`]
/// and the `with_*` setters:
///
/// ```
/// use openmassspec_core::{
///     Activation, Analyzer, Polarity, PrecursorInfo, ScanMode, SpectrumRecord,
/// };
///
/// let precursor = PrecursorInfo::default()
///     .with_target_mz(445.12)
///     .with_selected_mz(445.1201)
///     .with_isolation_width(1.6)
///     .with_charge(2)
///     .with_collision_energy(30.0)
///     .with_ce_is_nce(true)
///     .with_activation(Activation::HCD)
///     .with_precursor_native_id("controllerType=0 controllerNumber=1 scan=41".to_string());
///
/// let rec = SpectrumRecord::new(41, 42, "controllerType=0 controllerNumber=1 scan=42")
///     .with_ms_level(2)
///     .with_polarity(Polarity::Positive)
///     .with_scan_mode(ScanMode::Centroid)
///     .with_analyzer(Analyzer::FTMS)
///     .with_retention_time_sec(12.5 * 60.0)
///     .with_faims_cv(None) // `Option<T>` values are accepted too
///     .with_precursor(precursor)
///     .with_peaks(vec![150.1, 250.2], vec![10.0, 20.0])
///     .with_extra("myreader.segment", "2");
///
/// assert_eq!(rec.ms_level, 2);
/// assert_eq!(rec.precursor.as_ref().unwrap().charge, Some(2));
/// // Public fields stay readable and assignable.
/// let mut rec = rec;
/// rec.filter = Some("FTMS + c NSI d Full ms2 445.12@hcd30.00".into());
/// ```
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
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
    /// MS level (1 for MS1, 2 for MS2, ...).
    pub ms_level: u32,
    /// Scan polarity.
    pub polarity: Option<Polarity>,
    /// Centroid or profile.
    pub scan_mode: Option<ScanMode>,
    /// Mass analyzer that recorded this spectrum.
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
    /// Ion injection (accumulation) time in milliseconds.
    pub ion_injection_time_ms: Option<f64>,
    /// Mean inverse reduced ion mobility (1/K0) for the spectrum, when
    /// applicable (Bruker timsTOF, Waters TWIMS).
    pub inv_mobility: Option<f64>,
    /// FAIMS compensation voltage in volts, when the instrument has a FAIMS
    /// Pro interface and reports it per-scan (Thermo Orbitrap Fusion/Lumos/
    /// Eclipse/Ascend). `None` when not applicable or not decoded.
    pub faims_cv: Option<f64>,
    /// Precursor information for MS2+ spectra.
    pub precursor: Option<PrecursorInfo>,
    /// m/z values, ascending.
    pub mz: Vec<f64>,
    /// Intensities, parallel to `mz`.
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
    /// A record with the identifying fields set, `ms_level` 1, retention
    /// time 0, no peaks, and every optional field `None`. Fill in the rest
    /// with the `with_*` setters or by assigning the public fields.
    pub fn new(index: usize, scan_number: u32, native_id: impl Into<String>) -> Self {
        Self {
            index,
            scan_number,
            native_id: native_id.into(),
            ms_level: 1,
            polarity: None,
            scan_mode: None,
            analyzer: None,
            acquisition_event_id: None,
            filter: None,
            retention_time_sec: 0.0,
            total_ion_current: None,
            base_peak_mz: None,
            base_peak_intensity: None,
            low_mz: None,
            high_mz: None,
            ion_injection_time_ms: None,
            inv_mobility: None,
            faims_cv: None,
            precursor: None,
            mz: Vec::new(),
            intensity: Vec::new(),
            inv_mobility_per_peak: None,
            extra: BTreeMap::new(),
        }
    }

    /// Set the `mz` and `intensity` arrays together.
    pub fn with_peaks(mut self, mz: Vec<f64>, intensity: Vec<f32>) -> Self {
        self.mz = mz;
        self.intensity = intensity;
        self
    }

    /// Insert one `extra` entry. Keys use `<reader>.<field>`.
    pub fn with_extra(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.extra.insert(key.into(), value.into());
        self
    }

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
            _ => match (self.mz.first(), self.mz.last()) {
                (Some(&lo), Some(&hi)) => Some((lo, hi)),
                _ => None,
            },
        }
    }
}

/// A chromatogram trace (TIC, BPC, SRM/MRM transition).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct ChromatogramRecord {
    /// Zero-based position in the chromatogram list (mzML `index=`).
    pub index: usize,
    /// mzML chromatogram `id`.
    pub id: String,
    /// e.g. `"total ion current chromatogram"`, `"basepeak chromatogram"`,
    /// `"selected reaction monitoring chromatogram"`. The writer emits
    /// `MS:1000235` (total ion current chromatogram) when `None`.
    pub chromatogram_type: Option<CvTerm>,
    /// Precursor (Q1) m/z for SRM/MRM traces.
    pub precursor_mz: Option<f64>,
    /// Product (Q3) m/z for SRM/MRM traces.
    pub product_mz: Option<f64>,
    /// Retention time in seconds.
    pub time_sec: Vec<f32>,
    /// Intensities, parallel to `time_sec`.
    pub intensity: Vec<f32>,
}

impl ChromatogramRecord {
    /// A chromatogram with the given id and arrays and no type, precursor,
    /// or product.
    pub fn new(
        index: usize,
        id: impl Into<String>,
        time_sec: Vec<f32>,
        intensity: Vec<f32>,
    ) -> Self {
        Self {
            index,
            id: id.into(),
            chromatogram_type: None,
            precursor_mz: None,
            product_mz: None,
            time_sec,
            intensity,
        }
    }
}

/// Run-level metadata.
///
/// Vendors construct this once per file with [`RunMetadata::new`] and the
/// `with_*` setters. The mzML writer uses it to populate
/// `<fileDescription>`, `<sourceFileList>`, `<instrumentConfigurationList>`
/// and `<softwareList>`.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
#[non_exhaustive]
pub struct RunMetadata {
    /// The source file name to put in `<sourceFile name="...">`. Typically
    /// `Path::file_name()`; for directory-based formats (Bruker `.d/`,
    /// Waters `.raw/`), the directory name. The `<run id>` is derived from
    /// it (sanitized to a valid `xs:ID`).
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
    /// Software identifier (e.g. `"opentfraw"`). The `<software id>` is
    /// derived from it (sanitized to a valid `xs:ID`); the unmodified name
    /// is the value of the software's `MS:1000799` cvParam.
    pub software_name: String,
    /// Version of `software_name`.
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
    /// ([`MobilityArrayKind::InverseReducedVsPerCm2`]).
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
    /// a Fusion) preserve per-scan analyzer identity. When empty, every
    /// spectrum references the single default `IC1`.
    pub analyzers: Vec<Analyzer>,
    /// Namespaced reader-specific values. Keys use `<reader>.<field>`.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub extra: BTreeMap<String, String>,
}

impl RunMetadata {
    /// Metadata with the required fields set and every optional field
    /// `None` or empty.
    pub fn new(
        source_file_name: impl Into<String>,
        source_file_format: CvTerm,
        native_id_format: CvTerm,
        instrument: CvTerm,
        software_name: impl Into<String>,
        software_version: impl Into<String>,
    ) -> Self {
        Self {
            source_file_name: source_file_name.into(),
            source_file_format,
            native_id_format,
            instrument,
            instrument_serial_number: None,
            software_name: software_name.into(),
            software_version: software_version.into(),
            acquisition_software_name: None,
            acquisition_software_version: None,
            start_timestamp: None,
            mobility_array_kind: None,
            analyzers: Vec::new(),
            extra: BTreeMap::new(),
        }
    }

    /// Insert one `extra` entry. Keys use `<reader>.<field>`.
    pub fn with_extra(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.extra.insert(key.into(), value.into());
        self
    }
}

setters!(SpectrumRecord {
    val with_index index: usize,
    val with_scan_number scan_number: u32,
    into with_native_id native_id: String,
    val with_ms_level ms_level: u32,
    opt with_polarity polarity: Polarity,
    opt with_scan_mode scan_mode: ScanMode,
    opt with_analyzer analyzer: Analyzer,
    opt with_acquisition_event_id acquisition_event_id: u32,
    opt with_filter filter: String,
    val with_retention_time_sec retention_time_sec: f64,
    opt with_total_ion_current total_ion_current: f64,
    opt with_base_peak_mz base_peak_mz: f64,
    opt with_base_peak_intensity base_peak_intensity: f64,
    opt with_low_mz low_mz: f64,
    opt with_high_mz high_mz: f64,
    opt with_ion_injection_time_ms ion_injection_time_ms: f64,
    opt with_inv_mobility inv_mobility: f64,
    opt with_faims_cv faims_cv: f64,
    opt with_precursor precursor: PrecursorInfo,
    val with_mz mz: Vec<f64>,
    val with_intensity intensity: Vec<f32>,
    opt with_inv_mobility_per_peak inv_mobility_per_peak: Vec<f32>,
    val with_extra_map extra: BTreeMap<String, String>,
});

setters!(PrecursorInfo {
    opt with_target_mz target_mz: f64,
    opt with_selected_mz selected_mz: f64,
    opt with_isolation_width isolation_width: f64,
    opt with_charge charge: i32,
    opt with_intensity intensity: f64,
    opt with_collision_energy collision_energy: f64,
    val with_ce_is_nce ce_is_nce: bool,
    opt with_precursor_native_id precursor_native_id: String,
    opt with_activation activation: Activation,
    opt with_analyzer analyzer: Analyzer,
    opt with_ccs ccs: f64,
});

setters!(ChromatogramRecord {
    val with_index index: usize,
    into with_id id: String,
    opt with_chromatogram_type chromatogram_type: CvTerm,
    opt with_precursor_mz precursor_mz: f64,
    opt with_product_mz product_mz: f64,
    val with_time_sec time_sec: Vec<f32>,
    val with_intensity intensity: Vec<f32>,
});

setters!(RunMetadata {
    into with_source_file_name source_file_name: String,
    val with_source_file_format source_file_format: CvTerm,
    val with_native_id_format native_id_format: CvTerm,
    val with_instrument instrument: CvTerm,
    opt with_instrument_serial_number instrument_serial_number: String,
    into with_software_name software_name: String,
    into with_software_version software_version: String,
    opt with_acquisition_software_name acquisition_software_name: String,
    opt with_acquisition_software_version acquisition_software_version: String,
    opt with_start_timestamp start_timestamp: String,
    opt with_mobility_array_kind mobility_array_kind: MobilityArrayKind,
    val with_analyzers analyzers: Vec<Analyzer>,
    val with_extra_map extra: BTreeMap<String, String>,
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deserialized_cv_terms_own_their_strings() {
        let json = r#"{"accession":"MS:1000563","name":"Thermo RAW format"}"#;
        let a: CvTerm = serde_json::from_str(json).unwrap();
        assert_eq!(a, CvTerm::new("MS:1000563", "Thermo RAW format"));
        assert!(matches!(a.accession, Cow::Owned(_)));
        assert!(matches!(a.name, Cow::Owned(_)));
        let generated: CvTerm = cv::PEAK_PICKING.into();
        assert!(matches!(generated.accession, Cow::Borrowed("MS:1000035")));
    }

    fn spectrum() -> SpectrumRecord {
        SpectrumRecord::new(0, 1, "scan=1")
            .with_polarity(Polarity::Positive)
            .with_scan_mode(ScanMode::Centroid)
            .with_analyzer(Analyzer::TOFMS)
            .with_acquisition_event_id(17)
            .with_retention_time_sec(1.5)
            .with_peaks(vec![100.0], vec![12.0])
            .with_extra("openszraw.event_label", "survey")
    }

    fn run_metadata() -> RunMetadata {
        RunMetadata::new(
            "test.raw",
            CvTerm::new("MS:1000563", "Thermo RAW format"),
            CvTerm::new("MS:1000768", "Thermo nativeID format"),
            CvTerm::new("MS:1001911", "Q Exactive"),
            "test-reader",
            "1.0",
        )
        .with_extra("openwraw.source_type", "ESI")
    }

    #[test]
    fn setters_accept_values_and_options() {
        let rec = SpectrumRecord::new(3, 4, "scan=4")
            .with_faims_cv(-45.0)
            .with_low_mz(None)
            .with_filter(Some("f".to_string()));
        assert_eq!(rec.faims_cv, Some(-45.0));
        assert_eq!(rec.low_mz, None);
        assert_eq!(rec.filter.as_deref(), Some("f"));
        assert_eq!(rec.ms_level, 1);
    }

    #[test]
    fn serde_round_trips_records() {
        let s = spectrum();
        let json = serde_json::to_string(&s).unwrap();
        let back: SpectrumRecord = serde_json::from_str(&json).unwrap();
        assert_eq!(back, s);
        assert_eq!(back.extra["openszraw.event_label"], "survey");

        let m = run_metadata();
        let json = serde_json::to_string(&m).unwrap();
        let back: RunMetadata = serde_json::from_str(&json).unwrap();
        assert_eq!(back, m);

        let c = ChromatogramRecord::new(0, "TIC", vec![0.0, 1.0], vec![2.0, 3.0])
            .with_chromatogram_type(CvTerm::from(cv::TOTAL_ION_CURRENT_CHROMATOGRAM));
        let json = serde_json::to_string(&c).unwrap();
        let back: ChromatogramRecord = serde_json::from_str(&json).unwrap();
        assert_eq!(back, c);
    }

    #[test]
    fn serde_defaults_optional_fields_when_absent() {
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
