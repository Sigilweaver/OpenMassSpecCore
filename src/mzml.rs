//! Canonical mzML 1.1.0 writer.
//!
//! Consumes any [`SpectrumSource`] and emits valid mzML, with an indexed
//! variant that also writes the `<indexList>` + `<fileChecksum>` trailer so
//! random-access mzML readers (pyteomics, pymzml, sciex-style indexers) work.
//!
//! Binary data arrays (`m/z`, intensity, and per-peak ion mobility, plus
//! chromatogram time/intensity) are zlib-compressed (`MS:1000574`) by
//! default, matching msconvert/OpenMS output; set
//! [`WriteOptions::with_compression`] to write uncompressed (`MS:1000576`)
//! arrays instead.
//!
//! Record errors: by default the writer stops at the first `Err` the source
//! yields and returns it, leaving a truncated document in `out`. With
//! [`RecordErrorPolicy::Skip`] it leaves failed records out, counts them, and
//! reports the counts in the returned [`WriteSummary`].
//!
//! Every PSI-MS / UO term the writer emits on its own comes from the typed
//! constants in [`crate::cv`], so accession and name cannot disagree.
//! Vendor-specific terms (the source-file format, the native-ID format, the
//! instrument, chromatogram types) come from the source's
//! [`RunMetadata`](crate::RunMetadata) and records as [`CvTerm`] values.

use std::fmt::Display;
use std::io::{self, Write};

use flate2::write::ZlibEncoder;
use flate2::Compression as ZlibLevel;

use crate::cv::{self, Term};
use crate::enums::{Activation, Analyzer, Compression, MobilityArrayKind, Polarity, ScanMode};
use crate::error::Result;
use crate::source::SpectrumSource;
use crate::types::{ChromatogramRecord, CvTerm, RunMetadata, SpectrumRecord};

// ---------- byte-counting writer that also feeds a streaming SHA-1 ----------

struct CountingWriter<'a, W: Write> {
    inner: &'a mut W,
    pos: u64,
    sha1: Sha1,
    hashing: bool,
}

impl<'a, W: Write> CountingWriter<'a, W> {
    fn new(inner: &'a mut W) -> Self {
        Self {
            inner,
            pos: 0,
            sha1: Sha1::new(),
            hashing: true,
        }
    }
}

impl<W: Write> Write for CountingWriter<'_, W> {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        let n = self.inner.write(buf)?;
        self.pos += n as u64;
        if self.hashing {
            self.sha1.update(&buf[..n]);
        }
        Ok(n)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        self.inner.flush()
    }
}

// ---------- minimal SHA-1 (RFC 3174) so we don't pull in a crypto dep -------

struct Sha1 {
    state: [u32; 5],
    count: u64,
    buf: [u8; 64],
    buf_len: usize,
}

impl Sha1 {
    fn new() -> Self {
        Self {
            state: [0x67452301, 0xEFCDAB89, 0x98BADCFE, 0x10325476, 0xC3D2E1F0],
            count: 0,
            buf: [0u8; 64],
            buf_len: 0,
        }
    }

    fn update(&mut self, data: &[u8]) {
        let mut off = 0;
        while off < data.len() {
            let space = 64 - self.buf_len;
            let take = space.min(data.len() - off);
            self.buf[self.buf_len..self.buf_len + take].copy_from_slice(&data[off..off + take]);
            self.buf_len += take;
            self.count += take as u64;
            off += take;
            if self.buf_len == 64 {
                self.compress();
                self.buf_len = 0;
            }
        }
    }

    fn compress(&mut self) {
        let mut w = [0u32; 80];
        for (i, word) in w.iter_mut().enumerate().take(16) {
            *word = u32::from_be_bytes(self.buf[i * 4..i * 4 + 4].try_into().unwrap());
        }
        for i in 16..80 {
            w[i] = (w[i - 3] ^ w[i - 8] ^ w[i - 14] ^ w[i - 16]).rotate_left(1);
        }
        let [mut a, mut b, mut c, mut d, mut e] = self.state;
        for (i, &wi) in w.iter().enumerate() {
            let (f, k) = match i {
                0..=19 => ((b & c) | (!b & d), 0x5A827999u32),
                20..=39 => (b ^ c ^ d, 0x6ED9EBA1),
                40..=59 => ((b & c) | (b & d) | (c & d), 0x8F1BBCDC),
                _ => (b ^ c ^ d, 0xCA62C1D6),
            };
            let temp = a
                .rotate_left(5)
                .wrapping_add(f)
                .wrapping_add(e)
                .wrapping_add(k)
                .wrapping_add(wi);
            e = d;
            d = c;
            c = b.rotate_left(30);
            b = a;
            a = temp;
        }
        self.state[0] = self.state[0].wrapping_add(a);
        self.state[1] = self.state[1].wrapping_add(b);
        self.state[2] = self.state[2].wrapping_add(c);
        self.state[3] = self.state[3].wrapping_add(d);
        self.state[4] = self.state[4].wrapping_add(e);
    }

    fn finalize(mut self) -> [u8; 20] {
        let bit_count = self.count * 8;
        self.update(&[0x80]);
        while self.buf_len != 56 {
            self.update(&[0u8]);
        }
        self.update(&bit_count.to_be_bytes());
        let mut digest = [0u8; 20];
        for (i, &word) in self.state.iter().enumerate() {
            digest[i * 4..i * 4 + 4].copy_from_slice(&word.to_be_bytes());
        }
        digest
    }
}

// ---------- base64 (RFC 4648 sec 4, no wrapping) ----------------------------

const B64: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

fn base64_encode(data: &[u8]) -> String {
    let n = data.len();
    let mut out = Vec::with_capacity(n.div_ceil(3) * 4);
    let mut i = 0;
    while i + 2 < n {
        let b = ((data[i] as u32) << 16) | ((data[i + 1] as u32) << 8) | (data[i + 2] as u32);
        out.push(B64[((b >> 18) & 0x3f) as usize]);
        out.push(B64[((b >> 12) & 0x3f) as usize]);
        out.push(B64[((b >> 6) & 0x3f) as usize]);
        out.push(B64[(b & 0x3f) as usize]);
        i += 3;
    }
    if n - i == 2 {
        let b = ((data[i] as u32) << 16) | ((data[i + 1] as u32) << 8);
        out.push(B64[((b >> 18) & 0x3f) as usize]);
        out.push(B64[((b >> 12) & 0x3f) as usize]);
        out.push(B64[((b >> 6) & 0x3f) as usize]);
        out.push(b'=');
    } else if n - i == 1 {
        let b = (data[i] as u32) << 16;
        out.push(B64[((b >> 18) & 0x3f) as usize]);
        out.push(B64[((b >> 12) & 0x3f) as usize]);
        out.push(b'=');
        out.push(b'=');
    }
    String::from_utf8(out).expect("base64 output is ASCII")
}

// ---------- zlib (RFC 1950) + compression-aware array encoding ------------

/// Compress `data` as a zlib stream. In-memory `Vec<u8>` I/O never fails, so
/// the `expect`s here are infallible in practice.
fn zlib_compress(data: &[u8]) -> Vec<u8> {
    let mut enc = ZlibEncoder::new(Vec::new(), ZlibLevel::default());
    enc.write_all(data)
        .expect("in-memory zlib compression cannot fail");
    enc.finish()
        .expect("in-memory zlib compression cannot fail")
}

/// The `MS:1000574`/`MS:1000576` term for a `<binaryDataArray>`, matching
/// whatever bytes [`encode_bytes`] wrote.
fn compression_term(compression: Compression) -> Term {
    match compression {
        Compression::NoCompression => cv::NO_COMPRESSION,
        Compression::Zlib => cv::ZLIB_COMPRESSION,
    }
}

/// Base64-encode `data`, zlib-compressing it first when `compression`
/// requests it. Pair with [`compression_term`] so the emitted cvParam always
/// matches the bytes actually written.
fn encode_bytes(data: &[u8], compression: Compression) -> String {
    match compression {
        Compression::NoCompression => base64_encode(data),
        Compression::Zlib => base64_encode(&zlib_compress(data)),
    }
}

fn encode_f64_array(vals: &[f64], compression: Compression) -> String {
    let bytes: Vec<u8> = vals.iter().flat_map(|v| v.to_le_bytes()).collect();
    encode_bytes(&bytes, compression)
}

fn encode_f32_array(vals: &[f32], compression: Compression) -> String {
    let bytes: Vec<u8> = vals.iter().flat_map(|v| v.to_le_bytes()).collect();
    encode_bytes(&bytes, compression)
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

// ---------- cvParam emission -------------------------------------------------

/// `<cvParam>` for a generated term. `value` must already be XML-escaped.
fn cv_param<W: Write>(
    out: &mut W,
    indent: &str,
    term: Term,
    value: impl Display,
) -> io::Result<()> {
    writeln!(
        out,
        r#"{indent}<cvParam cvRef="{}" accession="{}" name="{}" value="{value}"/>"#,
        term.cv_ref(),
        term.accession(),
        term.name()
    )
}

/// `<cvParam>` for a generated term with a unit. `value` must already be
/// XML-escaped.
fn cv_param_with_unit<W: Write>(
    out: &mut W,
    indent: &str,
    term: Term,
    value: impl Display,
    unit: Term,
) -> io::Result<()> {
    writeln!(
        out,
        r#"{indent}<cvParam cvRef="{}" accession="{}" name="{}" value="{value}" unitCvRef="{}" unitAccession="{}" unitName="{}"/>"#,
        term.cv_ref(),
        term.accession(),
        term.name(),
        unit.cv_ref(),
        unit.accession(),
        unit.name()
    )
}

/// `<cvParam value="">` for a caller-supplied [`CvTerm`].
fn cv_term_param<W: Write>(out: &mut W, indent: &str, term: &CvTerm) -> io::Result<()> {
    let cv_ref = if term.accession.starts_with("UO:") {
        "UO"
    } else {
        "MS"
    };
    writeln!(
        out,
        r#"{indent}<cvParam cvRef="{cv_ref}" accession="{}" name="{}" value=""/>"#,
        escape(&term.accession),
        escape(&term.name)
    )
}

/// PSI-MS dissociation-method terms for an [`Activation`].
///
/// An empty slice means no PSI-MS term is a faithful match, so the writer
/// emits no dissociation-method cvParam rather than guessing.
fn activation_terms(act: Activation, analyzer: Option<Analyzer>) -> &'static [Term] {
    match act {
        Activation::HCD => &[cv::BEAM_TYPE_COLLISION_INDUCED_DISSOCIATION],
        Activation::ETD => &[cv::ELECTRON_TRANSFER_DISSOCIATION],
        // ETD plus the combined EThcD term, so readers that only look for
        // ETD keep working while the supplemental HCD step is recorded.
        Activation::EThcD => &[
            cv::ELECTRON_TRANSFER_DISSOCIATION,
            cv::ELECTRON_TRANSFER_HIGHER_ENERGY_COLLISION_DISSOCIATION,
        ],
        Activation::CID => match analyzer {
            Some(Analyzer::FTMS) => &[cv::BEAM_TYPE_COLLISION_INDUCED_DISSOCIATION],
            _ => &[cv::COLLISION_INDUCED_DISSOCIATION],
        },
        // MS:1000435 carries the exact synonyms "multiphoton dissociation"
        // and "MPD" in the PSI-MS CV.
        Activation::MPID => &[cv::PHOTODISSOCIATION],
        Activation::ECD => &[cv::ELECTRON_CAPTURE_DISSOCIATION],
        Activation::IRMPD => &[cv::INFRARED_MULTIPHOTON_DISSOCIATION],
        // "Proton transfer decay / activated-ion" has no matching PSI-MS
        // dissociation term, so emit nothing.
        Activation::PD => &[],
        Activation::PQD => &[cv::PULSED_Q_DISSOCIATION],
        Activation::UVPD => &[cv::ULTRAVIOLET_PHOTODISSOCIATION],
        Activation::SID => &[cv::SURFACE_INDUCED_DISSOCIATION],
    }
}

// ---------- options and summary ----------------------------------------------

/// What the writer does when the source yields an `Err` record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum RecordErrorPolicy {
    /// Stop and return the error. The document in `out` is truncated.
    #[default]
    Fail,
    /// Leave the failed record out and keep writing. Each skipped record is
    /// counted in [`WriteSummary`]. The `<spectrumList count>` attribute is
    /// written before iteration starts, from
    /// [`SpectrumSource::spectrum_count_hint`], so it still includes skipped
    /// spectra; compare it against [`WriteSummary::spectra_written`].
    Skip,
}

/// Options for [`write_mzml_with_options`] and
/// [`write_indexed_mzml_with_options`].
///
/// ```
/// use openmassspec_core::{Compression, RecordErrorPolicy, WriteOptions};
///
/// let opts = WriteOptions::new()
///     .with_compression(Compression::NoCompression)
///     .with_record_error_policy(RecordErrorPolicy::Skip);
/// assert_eq!(opts.record_error_policy, RecordErrorPolicy::Skip);
/// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct WriteOptions {
    /// Binary-data-array codec. Defaults to [`Compression::Zlib`].
    pub compression: Compression,
    /// Handling of `Err` records. Defaults to [`RecordErrorPolicy::Fail`].
    pub record_error_policy: RecordErrorPolicy,
}

impl WriteOptions {
    /// Default options: zlib compression, stop at the first record error.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the binary-data-array codec.
    #[must_use]
    pub fn with_compression(mut self, compression: Compression) -> Self {
        self.compression = compression;
        self
    }

    /// Set the record-error policy.
    #[must_use]
    pub fn with_record_error_policy(mut self, policy: RecordErrorPolicy) -> Self {
        self.record_error_policy = policy;
        self
    }
}

/// What a successful write produced.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct WriteSummary {
    /// `<spectrum>` elements written.
    pub spectra_written: usize,
    /// Spectra left out because the source yielded `Err` (only with
    /// [`RecordErrorPolicy::Skip`]).
    pub spectra_skipped: usize,
    /// `<chromatogram>` elements written.
    pub chromatograms_written: usize,
    /// Chromatograms left out because the source yielded `Err` (only with
    /// [`RecordErrorPolicy::Skip`]).
    pub chromatograms_skipped: usize,
}

impl WriteSummary {
    /// Total records skipped (spectra plus chromatograms).
    pub fn skipped(&self) -> usize {
        self.spectra_skipped + self.chromatograms_skipped
    }
}

// ---------- public entry points --------------------------------------------

/// Write the source's spectra as mzML 1.1.0 (un-indexed) with default
/// [`WriteOptions`]: zlib-compressed arrays, stop at the first record error.
///
/// # Errors
///
/// Returns the first `Err` the source yields, or an I/O error from `out`.
/// `out` then holds a truncated document.
pub fn write_mzml<S: SpectrumSource + ?Sized, W: Write>(
    src: &mut S,
    out: &mut W,
) -> Result<WriteSummary> {
    write_mzml_with_options(src, out, &WriteOptions::default())
}

/// Write the source's spectra as mzML 1.1.0 (un-indexed) with explicit
/// [`WriteOptions`].
///
/// # Errors
///
/// Returns an I/O error from `out`, or the first `Err` the source yields
/// unless `opts` selects [`RecordErrorPolicy::Skip`]. `out` then holds a
/// truncated document.
pub fn write_mzml_with_options<S: SpectrumSource + ?Sized, W: Write>(
    src: &mut S,
    out: &mut W,
    opts: &WriteOptions,
) -> Result<WriteSummary> {
    let mut cw = CountingWriter::new(out);
    cw.hashing = false;
    let body = write_body(src, &mut cw, opts, false)?;
    writeln!(cw, r#"</mzML>"#)?;
    Ok(body.summary)
}

/// Write the source's spectra as indexed mzML 1.1.0 (with `<indexList>` and
/// `<fileChecksum>` trailer) with default [`WriteOptions`]: zlib-compressed
/// arrays, stop at the first record error.
///
/// # Errors
///
/// Returns the first `Err` the source yields, or an I/O error from `out`.
/// `out` then holds a truncated document.
pub fn write_indexed_mzml<S: SpectrumSource + ?Sized, W: Write>(
    src: &mut S,
    out: &mut W,
) -> Result<WriteSummary> {
    write_indexed_mzml_with_options(src, out, &WriteOptions::default())
}

/// Write the source's spectra as indexed mzML 1.1.0 (with `<indexList>` and
/// `<fileChecksum>` trailer) with explicit [`WriteOptions`].
///
/// The `<fileChecksum>` SHA-1 is computed over the bytes actually written to
/// `out`, from the start of the document up to and including the
/// `<fileChecksum>` open tag, as mzML 1.1 specifies. Each `<offset>` and the
/// `<indexListOffset>` give the byte position of the element's `<`. Skipped
/// records have no index entry.
///
/// # Errors
///
/// Returns an I/O error from `out`, or the first `Err` the source yields
/// unless `opts` selects [`RecordErrorPolicy::Skip`]. `out` then holds a
/// truncated document.
pub fn write_indexed_mzml_with_options<S: SpectrumSource + ?Sized, W: Write>(
    src: &mut S,
    out: &mut W,
    opts: &WriteOptions,
) -> Result<WriteSummary> {
    let mut cw = CountingWriter::new(out);
    let body = write_body(src, &mut cw, opts, true)?;
    writeln!(cw, r#"  </mzML>"#)?;

    cw.write_all(b"  ")?;
    let index_list_offset = cw.pos;
    let n_indexes = 1 + usize::from(!body.chromatogram_offsets.is_empty());
    writeln!(cw, r#"<indexList count="{n_indexes}">"#)?;
    writeln!(cw, r#"    <index name="spectrum">"#)?;
    for (id, offset) in &body.spectrum_offsets {
        writeln!(
            cw,
            r#"      <offset idRef="{}">{}</offset>"#,
            escape(id),
            offset
        )?;
    }
    writeln!(cw, r#"    </index>"#)?;
    if !body.chromatogram_offsets.is_empty() {
        writeln!(cw, r#"    <index name="chromatogram">"#)?;
        for (id, offset) in &body.chromatogram_offsets {
            writeln!(
                cw,
                r#"      <offset idRef="{}">{}</offset>"#,
                escape(id),
                offset
            )?;
        }
        writeln!(cw, r#"    </index>"#)?;
    }
    writeln!(cw, r#"  </indexList>"#)?;

    writeln!(
        cw,
        r#"  <indexListOffset>{}</indexListOffset>"#,
        index_list_offset
    )?;
    // mzML 1.1: the SHA-1 covers the file from the start up to and
    // including the `<fileChecksum>` open tag (as pwiz and mzdata do).
    write!(cw, r#"  <fileChecksum>"#)?;
    cw.hashing = false;
    let digest = std::mem::replace(&mut cw.sha1, Sha1::new()).finalize();
    let hex: String = digest.iter().map(|b| format!("{:02x}", b)).collect();
    writeln!(cw, r#"{}</fileChecksum>"#, hex)?;
    writeln!(cw, r#"</indexedmzML>"#)?;
    Ok(body.summary)
}

/// Output of [`write_body`]: counts plus, for indexed output, the byte
/// offset of every written record.
struct Body {
    summary: WriteSummary,
    spectrum_offsets: Vec<(String, u64)>,
    chromatogram_offsets: Vec<(String, u64)>,
}

/// Write everything from the XML declaration through `</run>`.
fn write_body<S: SpectrumSource + ?Sized, W: Write>(
    src: &mut S,
    cw: &mut CountingWriter<'_, W>,
    opts: &WriteOptions,
    indexed: bool,
) -> Result<Body> {
    let meta = src.run_metadata();
    let ids = DocumentIds::new(&meta);
    let count = src.spectrum_count_hint().unwrap_or(0);
    let mobility_kind = meta.mobility_array_kind;
    let analyzer_ic_ids = analyzer_ic_table(&meta.analyzers);
    let extra_processing = src.additional_processing_steps();
    let compression = opts.compression;
    let skip = opts.record_error_policy == RecordErrorPolicy::Skip;

    write_prologue(cw, &meta, &ids, count, indexed, &extra_processing)?;

    let mut summary = WriteSummary::default();
    let mut spectrum_offsets = Vec::new();
    for item in src.iter_spectra() {
        let rec = match item {
            Ok(rec) => rec,
            Err(_) if skip => {
                summary.spectra_skipped += 1;
                continue;
            }
            Err(e) => return Err(e),
        };
        // Offsets point at the element's `<`, so write the indent first.
        cw.write_all(RECORD_INDENT.as_bytes())?;
        if indexed {
            spectrum_offsets.push((rec.native_id.clone(), cw.pos));
        }
        write_spectrum(cw, &rec, mobility_kind, &analyzer_ic_ids, compression)?;
        summary.spectra_written += 1;
    }
    writeln!(cw, r#"    </spectrumList>"#)?;

    let mut chroms = Vec::new();
    for item in src.iter_chromatograms() {
        match item {
            Ok(rec) => chroms.push(rec),
            Err(_) if skip => summary.chromatograms_skipped += 1,
            Err(e) => return Err(e),
        }
    }
    let mut chromatogram_offsets = Vec::new();
    if !chroms.is_empty() {
        writeln!(
            cw,
            r#"    <chromatogramList count="{}" defaultDataProcessingRef="dp1">"#,
            chroms.len()
        )?;
        for rec in &chroms {
            cw.write_all(RECORD_INDENT.as_bytes())?;
            if indexed {
                chromatogram_offsets.push((rec.id.clone(), cw.pos));
            }
            write_chromatogram(cw, rec, compression)?;
            summary.chromatograms_written += 1;
        }
        writeln!(cw, r#"    </chromatogramList>"#)?;
    }

    writeln!(cw, r#"  </run>"#)?;
    Ok(Body {
        summary,
        spectrum_offsets,
        chromatogram_offsets,
    })
}

// ---------- prologue / spectrum body ---------------------------------------

/// Indent in front of each `<spectrum>`/`<chromatogram>` open tag. Callers
/// write it before calling [`write_spectrum`]/[`write_chromatogram`] so the
/// indexed writer can record offsets that point at the element's `<`.
const RECORD_INDENT: &str = "      ";

/// Fixed `xs:ID` of the acquisition-software `<software>` entry.
const ACQUISITION_SOFTWARE_ID: &str = "acquisition_software";

/// The `xs:ID` values derived from run metadata.
struct DocumentIds {
    /// `<software id>` of the reader software, also every
    /// `processingMethod softwareRef`.
    software: String,
    /// `<run id>`.
    run: String,
}

impl DocumentIds {
    fn new(meta: &RunMetadata) -> Self {
        let software = unique_ncname(&meta.software_name, "software_", is_fixed_id);
        let run = unique_ncname(&meta.source_file_name, "run_", |id| {
            is_fixed_id(id) || id == software
        });
        Self { software, run }
    }
}

/// Ids the writer emits as fixed literals.
fn is_fixed_id(id: &str) -> bool {
    id == "sf1"
        || id == "dp1"
        || id == ACQUISITION_SOFTWARE_ID
        || (id.starts_with("IC") && id[2..].parse::<usize>().is_ok())
}

/// Turn `raw` into an XML NCName (as `xs:ID` requires) that `taken` rejects.
///
/// Characters outside `[A-Za-z0-9._-]` become `_`. `prefix` (which starts
/// with a letter) is prepended while the result does not start with a letter
/// or `_`, or is taken. Names that are already valid and free pass through
/// unchanged.
fn unique_ncname(raw: &str, prefix: &str, taken: impl Fn(&str) -> bool) -> String {
    let mut id: String = raw
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '.' | '-' | '_') {
                c
            } else {
                '_'
            }
        })
        .collect();
    loop {
        let starts_ok = id
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_alphabetic() || c == '_');
        if starts_ok && !taken(&id) {
            return id;
        }
        id.insert_str(0, prefix);
    }
}

fn write_prologue<W: Write>(
    out: &mut W,
    meta: &RunMetadata,
    ids: &DocumentIds,
    n_spectra: usize,
    indexed: bool,
    extra_processing: &[CvTerm],
) -> io::Result<()> {
    writeln!(out, r#"<?xml version="1.0" encoding="utf-8"?>"#)?;
    if indexed {
        writeln!(out, r#"<indexedmzML xmlns="http://psi.hupo.org/ms/mzml""#)?;
        writeln!(
            out,
            r#"             xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance""#
        )?;
        writeln!(
            out,
            r#"             xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd">"#
        )?;
        writeln!(
            out,
            r#"  <mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">"#
        )?;
    } else {
        writeln!(out, r#"<mzML xmlns="http://psi.hupo.org/ms/mzml""#)?;
        writeln!(
            out,
            r#"      xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance""#
        )?;
        // Non-indexed documents have a plain <mzML> root, so they reference
        // the base schema (matching version="1.1.0"), not the indexed
        // wrapper schema used for <indexedmzML> output.
        writeln!(
            out,
            r#"      xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd""#
        )?;
        writeln!(out, r#"      version="1.1.0">"#)?;
    }

    // The MS entry names the exact psi-ms.obo release `crate::cv` was
    // generated from; the UO terms the writer uses come from the same file.
    writeln!(out, r#"  <cvList count="2">"#)?;
    writeln!(
        out,
        r#"    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="{}" URI="{}"/>"#,
        cv::PSI_MS_VERSION,
        cv::PSI_MS_SOURCE_URL
    )?;
    writeln!(
        out,
        r#"    <cv id="UO" fullName="Unit Ontology" version="09:04:2014" URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>"#
    )?;
    writeln!(out, r#"  </cvList>"#)?;

    writeln!(out, r#"  <fileDescription>"#)?;
    writeln!(out, r#"    <fileContent>"#)?;
    cv_param(out, "      ", cv::MS1_SPECTRUM, "")?;
    cv_param(out, "      ", cv::MSN_SPECTRUM, "")?;
    writeln!(out, r#"    </fileContent>"#)?;
    writeln!(out, r#"    <sourceFileList count="1">"#)?;
    writeln!(
        out,
        r#"      <sourceFile id="sf1" name="{}" location="">"#,
        escape(&meta.source_file_name)
    )?;
    cv_term_param(out, "        ", &meta.source_file_format)?;
    cv_term_param(out, "        ", &meta.native_id_format)?;
    writeln!(out, r#"      </sourceFile>"#)?;
    writeln!(out, r#"    </sourceFileList>"#)?;
    writeln!(out, r#"  </fileDescription>"#)?;

    writeln!(
        out,
        r#"  <softwareList count="{}">"#,
        1 + usize::from(meta.acquisition_software_name.is_some())
    )?;
    writeln!(
        out,
        r#"    <software id="{}" version="{}">"#,
        ids.software,
        escape(&meta.software_version)
    )?;
    cv_param(
        out,
        "      ",
        cv::CUSTOM_UNRELEASED_SOFTWARE_TOOL,
        escape(&meta.software_name),
    )?;
    writeln!(out, r#"    </software>"#)?;
    if let Some(name) = &meta.acquisition_software_name {
        // A fixed id (rather than one derived from `name`): no
        // `processingMethod` references this software, so the id only has
        // to be unique and valid.
        writeln!(
            out,
            r#"    <software id="{ACQUISITION_SOFTWARE_ID}" version="{}">"#,
            escape(meta.acquisition_software_version.as_deref().unwrap_or(""))
        )?;
        // Vendor acquisition software each has its own PSI-MS term, but
        // mapping vendor-reported names to accessions would need every
        // string verified against the CV, so use the generic term and keep
        // the real name as the value.
        cv_param(
            out,
            "      ",
            cv::CUSTOM_UNRELEASED_SOFTWARE_TOOL,
            escape(name),
        )?;
        writeln!(out, r#"    </software>"#)?;
    }
    writeln!(out, r#"  </softwareList>"#)?;

    writeln!(
        out,
        r#"  <instrumentConfigurationList count="{}">"#,
        1 + meta.analyzers.len()
    )?;
    writeln!(out, r#"    <instrumentConfiguration id="IC1">"#)?;
    cv_term_param(out, "      ", &meta.instrument)?;
    if let Some(serial) = &meta.instrument_serial_number {
        cv_param(out, "      ", cv::INSTRUMENT_SERIAL_NUMBER, escape(serial))?;
    }
    writeln!(out, r#"    </instrumentConfiguration>"#)?;
    for (analyzer, id) in analyzer_ic_table(&meta.analyzers) {
        writeln!(out, r#"    <instrumentConfiguration id="{id}">"#)?;
        cv_term_param(out, "      ", &meta.instrument)?;
        cv_param(out, "      ", analyzer_term(analyzer), "")?;
        writeln!(out, r#"    </instrumentConfiguration>"#)?;
    }
    writeln!(out, r#"  </instrumentConfigurationList>"#)?;

    writeln!(out, r#"  <dataProcessingList count="1">"#)?;
    writeln!(out, r#"    <dataProcessing id="dp1">"#)?;
    writeln!(
        out,
        r#"      <processingMethod order="0" softwareRef="{}">"#,
        ids.software
    )?;
    cv_param(out, "        ", cv::CONVERSION_TO_MZML, "")?;
    writeln!(out, r#"      </processingMethod>"#)?;
    // Adapters that transform spectra before they reach the writer (e.g.
    // `Centroided` peak-picking) contribute additional steps here so the
    // file's processing history stays consistent with what was actually
    // done to the data (see `SpectrumSource::additional_processing_steps`).
    for (order, term) in extra_processing.iter().enumerate() {
        writeln!(
            out,
            r#"      <processingMethod order="{}" softwareRef="{}">"#,
            order + 1,
            ids.software
        )?;
        cv_term_param(out, "        ", term)?;
        writeln!(out, r#"      </processingMethod>"#)?;
    }
    writeln!(out, r#"    </dataProcessing>"#)?;
    writeln!(out, r#"  </dataProcessingList>"#)?;

    let run_id = &ids.run;
    match &meta.start_timestamp {
        Some(ts) => writeln!(
            out,
            r#"  <run id="{run_id}" defaultInstrumentConfigurationRef="IC1" defaultSourceFileRef="sf1" startTimeStamp="{}">"#,
            escape(ts)
        )?,
        None => writeln!(
            out,
            r#"  <run id="{run_id}" defaultInstrumentConfigurationRef="IC1" defaultSourceFileRef="sf1">"#
        )?,
    }
    for (name, value) in &meta.extra {
        writeln!(
            out,
            r#"    <userParam name="{}" value="{}"/>"#,
            escape(name),
            escape(value)
        )?;
    }
    writeln!(
        out,
        r#"    <spectrumList count="{}" defaultDataProcessingRef="dp1">"#,
        n_spectra
    )?;
    Ok(())
}

/// PSI-MS mass-analyzer-component term for each [`Analyzer`] variant.
///
/// `FTMS` covers both Orbitrap and FT-ICR instruments (see the enum's own
/// doc comment); this maps it to `orbitrap` since that is by far the more
/// common physical analyzer behind a "FTMS" scan-filter token across the
/// vendors this crate supports. FT-ICR instruments will be mislabeled.
fn analyzer_term(analyzer: Analyzer) -> Term {
    match analyzer {
        Analyzer::ITMS => cv::ION_TRAP,
        Analyzer::TQMS | Analyzer::SQMS => cv::QUADRUPOLE,
        Analyzer::TOFMS => cv::TIME_OF_FLIGHT,
        Analyzer::FTMS => cv::ORBITRAP,
        Analyzer::Sector => cv::MAGNETIC_SECTOR,
    }
}

/// Assigns `instrumentConfiguration` ids (`IC2`, `IC3`, ...) to
/// `meta.analyzers`, in order. `IC1` is reserved for the default
/// configuration built from `meta.instrument`.
fn analyzer_ic_table(analyzers: &[Analyzer]) -> Vec<(Analyzer, String)> {
    analyzers
        .iter()
        .enumerate()
        .map(|(i, a)| (*a, format!("IC{}", i + 2)))
        .collect()
}

fn write_spectrum<W: Write>(
    out: &mut W,
    rec: &SpectrumRecord,
    mobility_kind: Option<MobilityArrayKind>,
    analyzer_ic_ids: &[(Analyzer, String)],
    compression: Compression,
) -> io::Result<()> {
    const I: &str = "        ";
    let spectrum_type = if rec.ms_level <= 1 {
        cv::MS1_SPECTRUM
    } else {
        cv::MSN_SPECTRUM
    };
    let n_peaks = rec.mz.len();

    let ic_ref = rec.analyzer.and_then(|a| {
        analyzer_ic_ids
            .iter()
            .find(|(candidate, _)| *candidate == a)
            .map(|(_, id)| id.as_str())
    });
    writeln!(
        out,
        r#"<spectrum id="{sid}" index="{idx}" defaultArrayLength="{n}">"#,
        sid = escape(&rec.native_id),
        idx = rec.index,
        n = n_peaks
    )?;
    cv_param(out, I, cv::MS_LEVEL, rec.ms_level)?;
    cv_param(out, I, spectrum_type, "")?;

    match rec.scan_mode {
        Some(ScanMode::Centroid) => cv_param(out, I, cv::CENTROID_SPECTRUM, "")?,
        _ => cv_param(out, I, cv::PROFILE_SPECTRUM, "")?,
    }

    match rec.polarity {
        Some(Polarity::Positive) => cv_param(out, I, cv::POSITIVE_SCAN, "")?,
        Some(Polarity::Negative) => cv_param(out, I, cv::NEGATIVE_SCAN, "")?,
        None => {}
    }

    let tic = rec.effective_tic();
    let (bp_mz, bp_int) = rec.effective_base_peak().unwrap_or((0.0, 0.0));
    let (lo_mz, hi_mz) = rec.effective_mz_range().unwrap_or((0.0, 0.0));

    cv_param(out, I, cv::TOTAL_ION_CURRENT, format_args!("{tic:.6}"))?;
    cv_param(out, I, cv::BASE_PEAK_MZ, format_args!("{bp_mz:.6}"))?;
    cv_param(out, I, cv::BASE_PEAK_INTENSITY, format_args!("{bp_int:.6}"))?;
    cv_param(out, I, cv::LOWEST_OBSERVED_MZ, format_args!("{lo_mz:.6}"))?;
    cv_param(out, I, cv::HIGHEST_OBSERVED_MZ, format_args!("{hi_mz:.6}"))?;

    if let Some(event_id) = rec.acquisition_event_id {
        writeln!(
            out,
            r#"        <userParam name="acquisition event id" value="{event_id}"/>"#
        )?;
    }
    for (name, value) in &rec.extra {
        writeln!(
            out,
            r#"        <userParam name="{}" value="{}"/>"#,
            escape(name),
            escape(value)
        )?;
    }

    writeln!(out, r#"        <scanList count="1">"#)?;
    cv_param(out, "          ", cv::NO_COMBINATION, "")?;
    match ic_ref {
        Some(id) => writeln!(out, r#"          <scan instrumentConfigurationRef="{id}">"#)?,
        None => writeln!(out, r#"          <scan>"#)?,
    }

    const S: &str = "            ";
    if let Some(f) = rec.filter.as_deref() {
        if !f.is_empty() {
            cv_param(out, S, cv::FILTER_STRING, escape(f))?;
        }
    }

    // mzML stores RT in minutes by convention.
    let rt_min = rec.retention_time_sec / 60.0;
    cv_param_with_unit(
        out,
        S,
        cv::SCAN_START_TIME,
        format_args!("{rt_min:.6}"),
        cv::MINUTE,
    )?;

    if let Some(it) = rec.ion_injection_time_ms {
        cv_param_with_unit(
            out,
            S,
            cv::ION_INJECTION_TIME,
            format_args!("{it:.6}"),
            cv::MILLISECOND,
        )?;
    }

    if let Some(mob) = rec.inv_mobility {
        cv_param_with_unit(
            out,
            S,
            cv::INVERSE_REDUCED_ION_MOBILITY,
            format_args!("{mob:.6}"),
            cv::VOLT_SECOND_PER_SQUARE_CENTIMETER,
        )?;
    }

    if let Some(faims) = rec.faims_cv {
        cv_param_with_unit(
            out,
            S,
            cv::FAIMS_COMPENSATION_VOLTAGE,
            format_args!("{faims:.4}"),
            cv::VOLT,
        )?;
    }

    writeln!(out, r#"            <scanWindowList count="1">"#)?;
    writeln!(out, r#"              <scanWindow>"#)?;
    const W16: &str = "                ";
    cv_param_with_unit(
        out,
        W16,
        cv::SCAN_WINDOW_LOWER_LIMIT,
        format_args!("{lo_mz:.6}"),
        cv::MZ,
    )?;
    cv_param_with_unit(
        out,
        W16,
        cv::SCAN_WINDOW_UPPER_LIMIT,
        format_args!("{hi_mz:.6}"),
        cv::MZ,
    )?;
    writeln!(out, r#"              </scanWindow>"#)?;
    writeln!(out, r#"            </scanWindowList>"#)?;
    writeln!(out, r#"          </scan>"#)?;
    writeln!(out, r#"        </scanList>"#)?;

    if let Some(pre) = rec.precursor.as_ref() {
        writeln!(out, r#"        <precursorList count="1">"#)?;
        if let Some(ref nid) = pre.precursor_native_id {
            writeln!(
                out,
                r#"          <precursor spectrumRef="{}">"#,
                escape(nid)
            )?;
        } else {
            writeln!(out, r#"          <precursor>"#)?;
        }

        const P: &str = "              ";
        if pre.target_mz.is_some() || pre.isolation_width.is_some() {
            writeln!(out, r#"            <isolationWindow>"#)?;
            if let Some(mz) = pre.target_mz {
                cv_param_with_unit(
                    out,
                    P,
                    cv::ISOLATION_WINDOW_TARGET_MZ,
                    format_args!("{mz:.6}"),
                    cv::MZ,
                )?;
            }
            if let Some(w) = pre.isolation_width {
                let half = w / 2.0;
                cv_param_with_unit(
                    out,
                    P,
                    cv::ISOLATION_WINDOW_LOWER_OFFSET,
                    format_args!("{half:.6}"),
                    cv::MZ,
                )?;
                cv_param_with_unit(
                    out,
                    P,
                    cv::ISOLATION_WINDOW_UPPER_OFFSET,
                    format_args!("{half:.6}"),
                    cv::MZ,
                )?;
            }
            writeln!(out, r#"            </isolationWindow>"#)?;
        }

        if let Some(mz) = pre.selected_mz {
            writeln!(out, r#"            <selectedIonList count="1">"#)?;
            writeln!(out, r#"              <selectedIon>"#)?;
            cv_param_with_unit(
                out,
                W16,
                cv::SELECTED_ION_MZ,
                format_args!("{mz:.6}"),
                cv::MZ,
            )?;
            if let Some(z) = pre.charge {
                cv_param(out, W16, cv::CHARGE_STATE, z)?;
            }
            if let Some(i) = pre.intensity {
                cv_param(out, W16, cv::PEAK_INTENSITY, format_args!("{i:.6}"))?;
            }
            if let Some(ccs) = pre.ccs {
                cv_param_with_unit(
                    out,
                    W16,
                    cv::COLLISIONAL_CROSS_SECTIONAL_AREA,
                    format_args!("{ccs:.6}"),
                    cv::SQUARE_ANGSTROM,
                )?;
            }
            writeln!(out, r#"              </selectedIon>"#)?;
            writeln!(out, r#"            </selectedIonList>"#)?;
        }

        writeln!(out, r#"            <activation>"#)?;
        // Unknown activation: emit no dissociation method rather than
        // asserting CID, which would be wrong for ETD/UVPD/... sources.
        if let Some(act) = pre.activation {
            for term in activation_terms(act, pre.analyzer) {
                cv_param(out, P, *term, "")?;
            }
        }
        if let Some(e) = pre.collision_energy {
            if pre.ce_is_nce {
                cv_param_with_unit(
                    out,
                    P,
                    cv::NORMALIZED_COLLISION_ENERGY,
                    format_args!("{e:.2}"),
                    cv::PERCENT,
                )?;
            } else {
                cv_param_with_unit(
                    out,
                    P,
                    cv::COLLISION_ENERGY,
                    format_args!("{e:.2}"),
                    cv::ELECTRONVOLT,
                )?;
            }
        }
        writeln!(out, r#"            </activation>"#)?;
        writeln!(out, r#"          </precursor>"#)?;
        writeln!(out, r#"        </precursorList>"#)?;
    }

    if n_peaks > 0 {
        let compression_cv = compression_term(compression);
        let mz_b64 = encode_f64_array(&rec.mz, compression);
        let int_b64 = encode_f32_array(&rec.intensity, compression);
        let mobility_b64_opt = rec
            .inv_mobility_per_peak
            .as_ref()
            .filter(|v| v.len() == n_peaks)
            .map(|v| encode_f32_array(v, compression));
        let array_count = 2 + usize::from(mobility_b64_opt.is_some());

        writeln!(
            out,
            r#"        <binaryDataArrayList count="{array_count}">"#
        )?;
        write_binary_array(
            out,
            &mz_b64,
            cv::MZ_ARRAY,
            None,
            cv::FLOAT_64_BIT,
            compression_cv,
        )?;
        write_binary_array(
            out,
            &int_b64,
            cv::INTENSITY_ARRAY,
            None,
            cv::FLOAT_32_BIT,
            compression_cv,
        )?;
        if let Some(mobility_b64) = mobility_b64_opt {
            let (array, unit) = match mobility_kind {
                Some(MobilityArrayKind::DriftTimeMilliseconds) => {
                    (cv::RAW_ION_MOBILITY_ARRAY, cv::MILLISECOND)
                }
                Some(MobilityArrayKind::InverseReducedVsPerCm2) | None => (
                    cv::RAW_INVERSE_REDUCED_ION_MOBILITY_ARRAY,
                    cv::VOLT_SECOND_PER_SQUARE_CENTIMETER,
                ),
            };
            write_binary_array(
                out,
                &mobility_b64,
                array,
                Some(unit),
                cv::FLOAT_32_BIT,
                compression_cv,
            )?;
        }
        writeln!(out, r#"        </binaryDataArrayList>"#)?;
    }

    writeln!(out, r#"      </spectrum>"#)?;
    Ok(())
}

/// One `<binaryDataArray>`: the array-type term (with `unit` when given),
/// the precision term, the compression term, and the base64 payload.
fn write_binary_array<W: Write>(
    out: &mut W,
    b64: &str,
    array: Term,
    unit: Option<Term>,
    precision: Term,
    compression: Term,
) -> io::Result<()> {
    const B: &str = "            ";
    writeln!(
        out,
        r#"          <binaryDataArray encodedLength="{}">"#,
        b64.len()
    )?;
    match unit {
        Some(unit) => cv_param_with_unit(out, B, array, "", unit)?,
        None => cv_param(out, B, array, "")?,
    }
    cv_param(out, B, precision, "")?;
    cv_param(out, B, compression, "")?;
    writeln!(out, r#"            <binary>{b64}</binary>"#)?;
    writeln!(out, r#"          </binaryDataArray>"#)?;
    Ok(())
}

fn write_chromatogram<W: Write>(
    out: &mut W,
    rec: &ChromatogramRecord,
    compression: Compression,
) -> io::Result<()> {
    let n = rec.time_sec.len();

    writeln!(
        out,
        r#"<chromatogram id="{id}" index="{idx}" defaultArrayLength="{n}">"#,
        id = escape(&rec.id),
        idx = rec.index,
    )?;

    match &rec.chromatogram_type {
        Some(term) => cv_term_param(out, "        ", term)?,
        None => cv_param(out, "        ", cv::TOTAL_ION_CURRENT_CHROMATOGRAM, "")?,
    }

    if let Some(mz) = rec.precursor_mz {
        writeln!(out, r#"        <precursor>"#)?;
        writeln!(out, r#"          <isolationWindow>"#)?;
        cv_param_with_unit(
            out,
            "            ",
            cv::ISOLATION_WINDOW_TARGET_MZ,
            format_args!("{mz:.6}"),
            cv::MZ,
        )?;
        writeln!(out, r#"          </isolationWindow>"#)?;
        // `<activation>` is required inside `<precursor>`, but
        // `ChromatogramRecord` carries no activation method, so leave it
        // empty rather than asserting CID.
        writeln!(out, r#"          <activation>"#)?;
        writeln!(out, r#"          </activation>"#)?;
        writeln!(out, r#"        </precursor>"#)?;
    }

    if let Some(mz) = rec.product_mz {
        writeln!(out, r#"        <product>"#)?;
        writeln!(out, r#"          <isolationWindow>"#)?;
        cv_param_with_unit(
            out,
            "            ",
            cv::ISOLATION_WINDOW_TARGET_MZ,
            format_args!("{mz:.6}"),
            cv::MZ,
        )?;
        writeln!(out, r#"          </isolationWindow>"#)?;
        writeln!(out, r#"        </product>"#)?;
    }

    // binaryDataArrayList is required (not optional) on ChromatogramType,
    // unlike SpectrumType where it's dropped for zero-peak spectra.
    // mzML stores time arrays in minutes by convention, matching scan start
    // time; `time_sec` is stored in seconds per `ChromatogramRecord`'s docs.
    let time_min: Vec<f32> = rec.time_sec.iter().map(|&t| t / 60.0).collect();
    let compression_cv = compression_term(compression);
    let time_b64 = encode_f32_array(&time_min, compression);
    let int_b64 = encode_f32_array(&rec.intensity, compression);

    writeln!(out, r#"        <binaryDataArrayList count="2">"#)?;
    write_binary_array(
        out,
        &time_b64,
        cv::TIME_ARRAY,
        Some(cv::MINUTE),
        cv::FLOAT_32_BIT,
        compression_cv,
    )?;
    write_binary_array(
        out,
        &int_b64,
        cv::INTENSITY_ARRAY,
        None,
        cv::FLOAT_32_BIT,
        compression_cv,
    )?;
    writeln!(out, r#"        </binaryDataArrayList>"#)?;
    writeln!(out, r#"      </chromatogram>"#)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enums::Polarity;

    #[test]
    fn base64_rfc_vectors() {
        assert_eq!(base64_encode(b""), "");
        assert_eq!(base64_encode(b"f"), "Zg==");
        assert_eq!(base64_encode(b"fo"), "Zm8=");
        assert_eq!(base64_encode(b"foo"), "Zm9v");
        assert_eq!(base64_encode(b"foobar"), "Zm9vYmFy");
        assert_eq!(base64_encode(b"Man"), "TWFu");
    }

    struct ToySource {
        source_file_name: String,
        start_timestamp: Option<String>,
        instrument_serial_number: Option<String>,
        acquisition_software_name: Option<String>,
        acquisition_software_version: Option<String>,
        analyzers: Vec<Analyzer>,
        extra: std::collections::BTreeMap<String, String>,
        spectra: Vec<SpectrumRecord>,
        chroms: Vec<ChromatogramRecord>,
        extra_processing: Vec<CvTerm>,
        /// Positions in `spectra` that are yielded as `Err` instead.
        failing_spectra: Vec<usize>,
        /// Positions in `chroms` that are yielded as `Err` instead.
        failing_chroms: Vec<usize>,
    }

    fn decode_error(what: &str, i: usize) -> crate::Error {
        crate::Error::vendor(std::io::Error::other(format!("{what} {i} is corrupt")))
    }

    impl ToySource {
        fn new() -> Self {
            Self {
                source_file_name: "toy.raw".into(),
                start_timestamp: None,
                instrument_serial_number: None,
                acquisition_software_name: None,
                acquisition_software_version: None,
                analyzers: Vec::new(),
                extra: Default::default(),
                spectra: vec![minimal_spectrum(0, None)],
                chroms: Vec::new(),
                extra_processing: Vec::new(),
                failing_spectra: Vec::new(),
                failing_chroms: Vec::new(),
            }
        }
    }

    impl SpectrumSource for ToySource {
        fn run_metadata(&self) -> RunMetadata {
            RunMetadata {
                source_file_name: self.source_file_name.clone(),
                source_file_format: CvTerm::new("MS:1000563", "Thermo RAW format"),
                native_id_format: CvTerm::new("MS:1000768", "Thermo nativeID format"),
                instrument: CvTerm::new("MS:1001911", "Q Exactive"),
                instrument_serial_number: self.instrument_serial_number.clone(),
                acquisition_software_name: self.acquisition_software_name.clone(),
                acquisition_software_version: self.acquisition_software_version.clone(),
                software_name: "toy".into(),
                software_version: "0.0.0".into(),
                start_timestamp: self.start_timestamp.clone(),
                mobility_array_kind: None,
                analyzers: self.analyzers.clone(),
                extra: self.extra.clone(),
            }
        }

        fn iter_spectra<'a>(&'a mut self) -> crate::SpectrumIter<'a> {
            let failing = &self.failing_spectra;
            Box::new(self.spectra.iter().enumerate().map(move |(i, rec)| {
                if failing.contains(&i) {
                    Err(decode_error("spectrum", i))
                } else {
                    Ok(rec.clone())
                }
            }))
        }

        fn iter_chromatograms<'a>(&'a mut self) -> crate::ChromatogramIter<'a> {
            let failing = &self.failing_chroms;
            Box::new(self.chroms.iter().enumerate().map(move |(i, rec)| {
                if failing.contains(&i) {
                    Err(decode_error("chromatogram", i))
                } else {
                    Ok(rec.clone())
                }
            }))
        }

        fn spectrum_count_hint(&self) -> Option<usize> {
            Some(self.spectra.len())
        }

        fn additional_processing_steps(&self) -> Vec<CvTerm> {
            self.extra_processing.clone()
        }
    }

    fn minimal_spectrum(index: usize, faims_cv: Option<f64>) -> SpectrumRecord {
        SpectrumRecord {
            index,
            scan_number: (index + 1) as u32,
            native_id: format!("controllerType=0 controllerNumber=1 scan={}", index + 1),
            ms_level: 1,
            polarity: Some(Polarity::Positive),
            scan_mode: Some(ScanMode::Centroid),
            analyzer: None,
            acquisition_event_id: None,
            filter: None,
            retention_time_sec: index as f64,
            total_ion_current: None,
            base_peak_mz: None,
            base_peak_intensity: None,
            low_mz: None,
            high_mz: None,
            ion_injection_time_ms: None,
            inv_mobility: None,
            faims_cv,
            precursor: None,
            mz: vec![100.0],
            intensity: vec![1.0],
            inv_mobility_per_peak: None,
            extra: Default::default(),
        }
    }

    fn tic_chromatogram(index: usize) -> ChromatogramRecord {
        ChromatogramRecord {
            index,
            id: "TIC".into(),
            chromatogram_type: Some(CvTerm::new("MS:1000235", "total ion current chromatogram")),
            precursor_mz: None,
            product_mz: None,
            time_sec: vec![0.0, 60.0, 120.0],
            intensity: vec![10.0, 20.0, 15.0],
        }
    }

    fn srm_chromatogram(index: usize) -> ChromatogramRecord {
        ChromatogramRecord {
            index,
            id: "SRM SIC Q1=524.3 Q3=136.1".into(),
            chromatogram_type: Some(CvTerm::new(
                "MS:1001473",
                "selected reaction monitoring chromatogram",
            )),
            precursor_mz: Some(524.3),
            product_mz: Some(136.1),
            time_sec: vec![0.0, 30.0],
            intensity: vec![5.0, 8.0],
        }
    }

    #[test]
    fn start_timestamp_emitted_when_present() {
        let mut src = ToySource::new();
        src.start_timestamp = Some("2026-01-01T12:00:00Z".into());
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"startTimeStamp="2026-01-01T12:00:00Z""#));
    }

    #[test]
    fn start_timestamp_omitted_when_absent() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(!s.contains("startTimeStamp"));
    }

    #[test]
    fn instrument_serial_number_emitted_when_present() {
        let mut src = ToySource::new();
        src.instrument_serial_number = Some("SN12345".into());
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(
            r#"<cvParam cvRef="MS" accession="MS:1000529" name="instrument serial number" value="SN12345"/>"#
        ));
    }

    #[test]
    fn instrument_serial_number_omitted_when_absent() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(!s.contains("MS:1000529"));
    }

    #[test]
    fn acquisition_software_emitted_as_second_software_entry_when_present() {
        let mut src = ToySource::new();
        src.acquisition_software_name = Some("Xcalibur".into());
        src.acquisition_software_version = Some("4.4.16.14".into());
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"<softwareList count="2">"#));
        assert!(s.contains(r#"<software id="acquisition_software" version="4.4.16.14">"#));
        assert!(s.contains(
            r#"<cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="Xcalibur"/>"#
        ));
        // Also on the indexed path.
        let mut buf = Vec::new();
        write_indexed_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"<software id="acquisition_software" version="4.4.16.14">"#));
    }

    #[test]
    fn acquisition_software_omitted_when_absent() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"<softwareList count="1">"#));
        assert!(!s.contains("acquisition_software"));
    }

    #[test]
    fn acquisition_software_version_defaults_to_empty_string_when_absent() {
        let mut src = ToySource::new();
        src.acquisition_software_name = Some("MassHunter".into());
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"<software id="acquisition_software" version="">"#));
    }

    fn is_ncname(s: &str) -> bool {
        let mut chars = s.chars();
        chars
            .next()
            .is_some_and(|c| c.is_ascii_alphabetic() || c == '_')
            && chars.all(|c| c.is_ascii_alphanumeric() || matches!(c, '.' | '-' | '_'))
    }

    fn run_id(meta: &RunMetadata) -> String {
        DocumentIds::new(meta).run
    }

    fn software_id(name: &str) -> String {
        let meta = RunMetadata {
            software_name: name.into(),
            ..ToySource::new().run_metadata()
        };
        DocumentIds::new(&meta).software
    }

    #[test]
    fn software_id_is_a_valid_ncname() {
        assert_eq!(software_id("opentfraw"), "opentfraw");
        assert_eq!(software_id("open tf raw"), "open_tf_raw");
        assert_eq!(software_id("2raw"), "software_2raw");
        assert_eq!(software_id(""), "software_");
        assert_eq!(software_id("sf1"), "software_sf1");
        assert_eq!(software_id("IC1"), "software_IC1");
        assert_eq!(
            software_id("acquisition_software"),
            "software_acquisition_software"
        );
        for name in ["0", "a b", "r\u{e9}ader", "a:b", "x&y", ""] {
            let id = software_id(name);
            assert!(is_ncname(&id), "{name:?} -> {id:?}");
        }
    }

    #[test]
    fn software_id_and_refs_use_sanitized_name_and_value_keeps_original() {
        struct Named(ToySource);
        impl SpectrumSource for Named {
            fn run_metadata(&self) -> RunMetadata {
                RunMetadata {
                    software_name: "1 reader & co".into(),
                    ..self.0.run_metadata()
                }
            }
            fn iter_spectra<'a>(&'a mut self) -> crate::SpectrumIter<'a> {
                self.0.iter_spectra()
            }
        }
        let mut src = Named(ToySource::new());
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"<software id="software_1_reader___co" version="0.0.0">"#));
        assert!(s.contains(r#"softwareRef="software_1_reader___co""#));
        assert!(s.contains(r#"name="custom unreleased software tool" value="1 reader &amp; co""#));
    }

    #[test]
    fn run_id_does_not_collide_with_software_id() {
        let meta = RunMetadata {
            source_file_name: "reader".into(),
            software_name: "reader".into(),
            ..ToySource::new().run_metadata()
        };
        let ids = DocumentIds::new(&meta);
        assert_eq!(ids.software, "reader");
        assert_eq!(ids.run, "run_reader");
    }

    #[test]
    fn run_id_is_a_valid_ncname() {
        let meta = |name: &str| RunMetadata {
            source_file_name: name.into(),
            ..ToySource::new().run_metadata()
        };
        assert_eq!(run_id(&meta("toy.raw")), "toy.raw");
        assert_eq!(run_id(&meta("_x.d")), "_x.d");
        assert_eq!(
            run_id(&meta("20260101_sample.raw")),
            "run_20260101_sample.raw"
        );
        assert_eq!(run_id(&meta("my run (1).wiff")), "my_run__1_.wiff");
        assert_eq!(run_id(&meta("-a")), "run_-a");
        assert_eq!(run_id(&meta("")), "run_");
        assert_eq!(run_id(&meta("dp1")), "run_dp1");
        assert_eq!(run_id(&meta("IC2")), "run_IC2");
        assert_eq!(run_id(&meta("toy")), "run_toy");
        for name in ["0.raw", "a b", "\u{e9}t\u{e9}.raw", "a:b", ""] {
            let id = run_id(&meta(name));
            assert!(is_ncname(&id), "{name:?} -> {id:?}");
        }
    }

    #[test]
    fn run_id_uses_sanitized_name_for_digit_leading_file() {
        let mut src = ToySource::new();
        src.source_file_name = "2026 sample.raw".into();
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"<run id="run_2026_sample.raw" "#));
        // sourceFile name keeps the original, unsanitized file name.
        assert!(s.contains(r#"<sourceFile id="sf1" name="2026 sample.raw""#));
    }

    #[test]
    fn no_extra_instrument_configurations_when_analyzers_empty() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"<instrumentConfigurationList count="1">"#));
        assert!(!s.contains("instrumentConfigurationRef"));
    }

    #[test]
    fn no_additional_processing_steps_leaves_data_processing_list_untouched() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"<dataProcessingList count="1">"#));
        assert!(s.contains(r#"<processingMethod order="0" softwareRef="toy">"#));
        assert!(!s.contains(r#"<processingMethod order="1""#));
    }

    #[test]
    fn additional_processing_steps_appended_in_order_after_conversion_step() {
        let mut src = ToySource::new();
        src.extra_processing = vec![
            cv::PEAK_PICKING.into(),
            CvTerm::new("MS:1000592", "smoothing"),
        ];
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();

        assert!(s.contains(r#"<dataProcessingList count="1">"#));
        let conversion = s.find(r#"<processingMethod order="0""#).unwrap();
        let peak_picking = s.find(r#"<processingMethod order="1""#).unwrap();
        let smoothing = s.find(r#"<processingMethod order="2""#).unwrap();
        assert!(conversion < peak_picking);
        assert!(peak_picking < smoothing);
        assert!(s.contains(
            r#"<cvParam cvRef="MS" accession="MS:1000035" name="peak picking" value=""/>"#
        ));
        assert!(
            s.contains(r#"<cvParam cvRef="MS" accession="MS:1000592" name="smoothing" value=""/>"#)
        );
    }

    #[test]
    fn per_spectrum_instrument_configuration_ref_matches_analyzer() {
        let mut src = ToySource::new();
        src.analyzers = vec![Analyzer::FTMS, Analyzer::ITMS];
        src.spectra = vec![
            SpectrumRecord {
                analyzer: Some(Analyzer::FTMS),
                ..minimal_spectrum(0, None)
            },
            SpectrumRecord {
                index: 1,
                native_id: "controllerType=0 controllerNumber=1 scan=2".into(),
                analyzer: Some(Analyzer::ITMS),
                ..minimal_spectrum(1, None)
            },
        ];
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();

        assert!(s.contains(r#"<instrumentConfigurationList count="3">"#));
        assert!(s.contains(r#"<instrumentConfiguration id="IC2">"#));
        assert!(
            s.contains(r#"<cvParam cvRef="MS" accession="MS:1000484" name="orbitrap" value=""/>"#)
        );
        assert!(s.contains(r#"<instrumentConfiguration id="IC3">"#));
        assert!(
            s.contains(r#"<cvParam cvRef="MS" accession="MS:1000264" name="ion trap" value=""/>"#)
        );

        assert!(s.contains(r#"instrumentConfigurationRef="IC2""#));
        assert!(s.contains(r#"instrumentConfigurationRef="IC3""#));
    }

    #[test]
    fn spectrum_with_unmatched_analyzer_gets_no_ref() {
        let mut src = ToySource::new();
        src.analyzers = vec![Analyzer::FTMS];
        src.spectra = vec![SpectrumRecord {
            analyzer: Some(Analyzer::TOFMS),
            ..minimal_spectrum(0, None)
        }];
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(!s.contains("instrumentConfigurationRef"));
    }

    #[test]
    fn faims_cv_emitted_when_present() {
        let mut src = ToySource::new();
        src.spectra = vec![minimal_spectrum(0, Some(-45.0))];
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(
            r#"accession="MS:1001581" name="FAIMS compensation voltage" value="-45.0000""#
        ));
    }

    #[test]
    fn faims_cv_omitted_when_absent() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(!s.contains("MS:1001581"));
    }

    #[test]
    fn ccs_emitted_on_selected_ion_when_present() {
        let mut src = ToySource::new();
        src.spectra = vec![SpectrumRecord {
            ms_level: 2,
            precursor: Some(crate::types::PrecursorInfo {
                selected_mz: Some(500.5),
                ccs: Some(153.4),
                ..Default::default()
            }),
            ..minimal_spectrum(0, None)
        }];
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(
            r#"<cvParam cvRef="MS" accession="MS:1002954" name="collisional cross sectional area" value="153.400000" unitCvRef="UO" unitAccession="UO:0000324" unitName="square angstrom"/>"#
        ));
    }

    #[test]
    fn ccs_omitted_when_absent() {
        let mut src = ToySource::new();
        src.spectra = vec![SpectrumRecord {
            ms_level: 2,
            precursor: Some(crate::types::PrecursorInfo {
                selected_mz: Some(500.5),
                ..Default::default()
            }),
            ..minimal_spectrum(0, None)
        }];
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(!s.contains("MS:1002954"));
    }

    fn ms2_with(activation: Option<Activation>, nce: Option<f64>) -> String {
        let mut src = ToySource::new();
        src.spectra = vec![SpectrumRecord {
            ms_level: 2,
            precursor: Some(crate::types::PrecursorInfo {
                selected_mz: Some(500.5),
                activation,
                collision_energy: nce,
                ce_is_nce: nce.is_some(),
                ..Default::default()
            }),
            ..minimal_spectrum(0, None)
        }];
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        String::from_utf8(buf).unwrap()
    }

    fn activation_block(xml: &str) -> &str {
        let start = xml.find("<activation>").unwrap();
        let end = xml.find("</activation>").unwrap();
        &xml[start..end]
    }

    #[test]
    fn normalized_collision_energy_uses_ms_1000138_with_percent_unit() {
        let s = ms2_with(Some(Activation::HCD), Some(30.0));
        assert!(s.contains(
            r#"<cvParam cvRef="MS" accession="MS:1000138" name="normalized collision energy" value="30.00" unitCvRef="UO" unitAccession="UO:0000187" unitName="percent"/>"#
        ));
        assert!(!s.contains("MS:1002013"));
    }

    #[test]
    fn unknown_activation_emits_no_dissociation_method() {
        let s = ms2_with(None, None);
        let block = activation_block(&s);
        assert!(!block.contains("cvParam"), "{block}");
    }

    #[test]
    fn activation_terms_match_psi_ms() {
        let cases: &[(Activation, &[&str])] = &[
            (Activation::SID, &["MS:1000136"]),
            (Activation::MPID, &["MS:1000435"]),
            (Activation::EThcD, &["MS:1000598", "MS:1002631"]),
            (Activation::PD, &[]),
        ];
        for (act, expected) in cases {
            let s = ms2_with(Some(*act), None);
            let block = activation_block(&s);
            assert_eq!(block.matches("<cvParam").count(), expected.len(), "{act:?}");
            for acc in *expected {
                assert!(block.contains(acc), "{act:?} missing {acc}");
            }
        }
    }

    #[test]
    fn chromatogram_list_omitted_when_source_yields_none() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(!s.contains("chromatogramList"));
        assert!(!s.contains("<chromatogram "));
    }

    #[test]
    fn chromatogram_list_emitted_with_tic_and_srm() {
        let mut src = ToySource::new();
        src.chroms = vec![tic_chromatogram(0), srm_chromatogram(1)];
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"<chromatogramList count="2" defaultDataProcessingRef="dp1">"#));
        assert!(s.contains(r#"<chromatogram id="TIC" index="0" defaultArrayLength="3">"#));
        assert!(s.contains(
            r#"<chromatogram id="SRM SIC Q1=524.3 Q3=136.1" index="1" defaultArrayLength="2">"#
        ));
        assert!(s.contains(r#"accession="MS:1000235" name="total ion current chromatogram""#));
        assert!(s.contains(
            r#"accession="MS:1001473" name="selected reaction monitoring chromatogram""#
        ));
        // SRM carries precursor/product isolation windows; TIC carries neither.
        assert!(s.contains("<precursor>"));
        assert!(s.contains("<product>"));
        // chromatogramList must close before </run>, after </spectrumList>.
        let spectrum_list_end = s.find("</spectrumList>").unwrap();
        let chrom_list_start = s.find("<chromatogramList").unwrap();
        let run_end = s.find("</run>").unwrap();
        assert!(spectrum_list_end < chrom_list_start);
        assert!(chrom_list_start < run_end);
    }

    #[test]
    fn indexed_mzml_adds_chromatogram_index_block() {
        let mut src = ToySource::new();
        src.chroms = vec![tic_chromatogram(0)];
        let mut buf = Vec::new();
        write_indexed_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"<indexList count="2">"#));
        assert!(s.contains(r#"<index name="spectrum">"#));
        assert!(s.contains(r#"<index name="chromatogram">"#));
        assert!(s.contains(r#"idRef="TIC""#));
    }

    #[test]
    fn indexed_mzml_has_single_index_block_without_chromatograms() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_indexed_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"<indexList count="1">"#));
        assert!(!s.contains(r#"<index name="chromatogram">"#));
    }

    #[test]
    fn index_offsets_point_at_element_open_tags() {
        let mut src = ToySource::new();
        src.spectra.push(SpectrumRecord {
            native_id: "controllerType=0 controllerNumber=1 scan=2".into(),
            ..minimal_spectrum(1, None)
        });
        src.chroms = vec![tic_chromatogram(0), srm_chromatogram(1)];
        let mut buf = Vec::new();
        write_indexed_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();

        let mut checked = 0;
        for line in s.lines() {
            let Some(rest) = line.trim().strip_prefix("<offset idRef=\"") else {
                continue;
            };
            let (id, rest) = rest.split_once('"').unwrap();
            let offset: usize = rest
                .trim_start_matches('>')
                .trim_end_matches("</offset>")
                .parse()
                .unwrap();
            let tail = &s[offset..];
            assert!(
                tail.starts_with(&format!("<spectrum id=\"{id}\""))
                    || tail.starts_with(&format!("<chromatogram id=\"{id}\"")),
                "offset {offset} for {id:?} points at {:?}",
                &tail[..20.min(tail.len())]
            );
            checked += 1;
        }
        assert_eq!(checked, 4);

        let start = s.find("<indexListOffset>").unwrap() + "<indexListOffset>".len();
        let end = s.find("</indexListOffset>").unwrap();
        let offset: usize = s[start..end].parse().unwrap();
        assert!(s[offset..].starts_with("<indexList "));
    }

    // ---------- record errors ------------------------------------------

    #[test]
    fn record_error_fails_the_write_by_default() {
        let mut src = ToySource::new();
        src.spectra.push(minimal_spectrum(1, None));
        src.failing_spectra = vec![1];
        let mut buf = Vec::new();
        let err = write_mzml(&mut src, &mut buf).unwrap_err();
        assert!(matches!(err, crate::Error::Vendor(_)), "{err:?}");
        assert!(err.to_string().contains("spectrum 1 is corrupt"));

        let mut buf = Vec::new();
        let err = write_indexed_mzml(&mut src, &mut buf).unwrap_err();
        assert!(err.to_string().contains("spectrum 1 is corrupt"));
        // The document stops at the failing record: no trailer.
        assert!(!String::from_utf8(buf).unwrap().contains("<indexList"));
    }

    #[test]
    fn chromatogram_error_fails_the_write_by_default() {
        let mut src = ToySource::new();
        src.chroms = vec![tic_chromatogram(0)];
        src.failing_chroms = vec![0];
        let err = write_mzml(&mut src, &mut Vec::new()).unwrap_err();
        assert!(err.to_string().contains("chromatogram 0 is corrupt"));
    }

    #[test]
    fn skip_policy_counts_and_reports_skipped_records() {
        let mut src = ToySource::new();
        src.spectra = (0..4).map(|i| minimal_spectrum(i, None)).collect();
        src.failing_spectra = vec![1, 3];
        src.chroms = vec![tic_chromatogram(0), srm_chromatogram(1)];
        src.failing_chroms = vec![0];
        let opts = WriteOptions::new().with_record_error_policy(RecordErrorPolicy::Skip);

        for indexed in [false, true] {
            let mut buf = Vec::new();
            let summary = if indexed {
                write_indexed_mzml_with_options(&mut src, &mut buf, &opts).unwrap()
            } else {
                write_mzml_with_options(&mut src, &mut buf, &opts).unwrap()
            };
            assert_eq!(summary.spectra_written, 2);
            assert_eq!(summary.spectra_skipped, 2);
            assert_eq!(summary.chromatograms_written, 1);
            assert_eq!(summary.chromatograms_skipped, 1);
            assert_eq!(summary.skipped(), 3);

            let s = String::from_utf8(buf).unwrap();
            assert_eq!(s.matches("<spectrum ").count(), 2);
            assert!(s.contains(r#"<chromatogramList count="1""#));
            assert!(!s.contains(r#"id="TIC""#));
            if indexed {
                assert_eq!(s.matches("<offset ").count(), 3);
            }
        }
    }

    #[test]
    fn summary_counts_everything_when_nothing_fails() {
        let mut src = ToySource::new();
        src.chroms = vec![tic_chromatogram(0)];
        let summary = write_indexed_mzml(&mut src, &mut Vec::new()).unwrap();
        assert_eq!(summary.spectra_written, 1);
        assert_eq!(summary.chromatograms_written, 1);
        assert_eq!(summary.skipped(), 0);
    }

    #[test]
    fn generated_term_names_need_no_xml_escaping() {
        for t in cv::ALL {
            assert_eq!(escape(t.name()), t.name(), "{}", t.accession());
        }
    }

    #[test]
    fn cv_list_names_the_generated_obo_release() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(&format!(
            r#"version="{}" URI="{}""#,
            cv::PSI_MS_VERSION,
            cv::PSI_MS_SOURCE_URL
        )));
    }

    // ---------- compression ---------------------------------------------

    #[test]
    fn compression_term_selects_matching_accession() {
        assert_eq!(
            compression_term(Compression::NoCompression),
            cv::NO_COMPRESSION
        );
        assert_eq!(compression_term(Compression::Zlib), cv::ZLIB_COMPRESSION);
    }

    #[test]
    fn zlib_compress_round_trips() {
        use std::io::Read;
        let original = b"the quick brown fox jumps over the lazy dog, repeated: \
                          the quick brown fox jumps over the lazy dog";
        let compressed = zlib_compress(original);
        assert!(
            compressed.len() < original.len(),
            "zlib should shrink a repetitive payload"
        );
        let mut decoder = flate2::read::ZlibDecoder::new(&compressed[..]);
        let mut round_tripped = Vec::new();
        decoder.read_to_end(&mut round_tripped).unwrap();
        assert_eq!(round_tripped, original);
    }

    #[test]
    fn default_writer_emits_zlib_compression_not_no_compression() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"accession="MS:1000574" name="zlib compression""#));
        assert!(!s.contains("MS:1000576"));
    }

    #[test]
    fn indexed_default_writer_also_emits_zlib_compression() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_indexed_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"accession="MS:1000574" name="zlib compression""#));
        assert!(!s.contains("MS:1000576"));
    }

    #[test]
    fn explicit_no_compression_writes_raw_little_endian_arrays() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        let opts = WriteOptions::new().with_compression(Compression::NoCompression);
        write_mzml_with_options(&mut src, &mut buf, &opts).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(r#"accession="MS:1000576" name="no compression""#));
        assert!(!s.contains("MS:1000574"));
        // The single toy spectrum has mz = [100.0]; with no compression the
        // <binary> payload is exactly the base64 of its 8 little-endian
        // bytes.
        let expected = base64_encode(&100.0f64.to_le_bytes());
        assert!(s.contains(&format!("<binary>{expected}</binary>")));
    }

    /// Extracts the text between the first `<binary>...</binary>` pair in
    /// `xml`, matching the m/z array (always the first `binaryDataArray` the
    /// writer emits for a non-empty spectrum).
    fn first_binary_payload(xml: &str) -> &str {
        let start = xml.find("<binary>").expect("no <binary> element") + "<binary>".len();
        let end = xml[start..]
            .find("</binary>")
            .expect("unterminated <binary>")
            + start;
        &xml[start..end]
    }

    /// Minimal RFC 4648 base64 decoder, the inverse of this module's
    /// `base64_encode`, scoped to tests so round-trip assertions don't need
    /// an extra dependency.
    fn base64_decode(s: &str) -> Vec<u8> {
        fn val(c: u8) -> u8 {
            match c {
                b'A'..=b'Z' => c - b'A',
                b'a'..=b'z' => c - b'a' + 26,
                b'0'..=b'9' => c - b'0' + 52,
                b'+' => 62,
                b'/' => 63,
                _ => 0,
            }
        }
        let bytes: Vec<u8> = s.bytes().filter(|&b| b != b'=').collect();
        let mut out = Vec::new();
        let mut i = 0;
        while i < bytes.len() {
            let n = (bytes.len() - i).min(4);
            let mut chunk = [0u8; 4];
            chunk[..n].copy_from_slice(&bytes[i..i + n]);
            let b = ((val(chunk[0]) as u32) << 18)
                | ((val(chunk[1]) as u32) << 12)
                | ((val(chunk[2]) as u32) << 6)
                | (val(chunk[3]) as u32);
            if n > 1 {
                out.push((b >> 16) as u8);
            }
            if n > 2 {
                out.push((b >> 8) as u8);
            }
            if n > 3 {
                out.push(b as u8);
            }
            i += 4;
        }
        out
    }

    #[test]
    fn zlib_binary_payload_decodes_back_to_original_mz_array() {
        use std::io::Read;
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();

        let payload = first_binary_payload(&s);
        let compressed = base64_decode(payload);
        let mut decoder = flate2::read::ZlibDecoder::new(&compressed[..]);
        let mut raw = Vec::new();
        decoder.read_to_end(&mut raw).unwrap();

        // ToySource's single spectrum carries mz = [100.0] (see
        // `minimal_spectrum`), encoded little-endian as f64.
        assert_eq!(raw, 100.0f64.to_le_bytes());
    }

    #[test]
    fn indexed_zlib_checksum_still_present_and_well_formed() {
        let mut src = ToySource::new();
        let mut buf = Vec::new();
        write_indexed_mzml(&mut src, &mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // fileChecksum is hashed over whatever bytes actually got written
        // (see write_indexed_mzml_with_compression's doc comment), so the
        // trailer shape doesn't change with the codec - only the byte
        // count feeding the hash does.
        assert!(s.contains("<indexListOffset>"));
        let checksum_start = s.find("<fileChecksum>").unwrap() + "<fileChecksum>".len();
        let checksum_end = s.find("</fileChecksum>").unwrap();
        let hex = &s[checksum_start..checksum_end];
        assert_eq!(hex.len(), 40, "SHA-1 hex digest must be 40 chars");
        assert!(hex.chars().all(|c| c.is_ascii_hexdigit()));
    }

    #[test]
    fn projects_vendor_extras_and_acquisition_event_id_as_user_params() {
        let mut src = ToySource::new();
        src.extra
            .insert("openwraw.source_type".into(), "ESI & APCI".into());
        src.spectra[0].acquisition_event_id = Some(42);
        src.spectra[0]
            .extra
            .insert("openszraw.segment".into(), "A < B".into());

        let mut buf = Vec::new();
        write_mzml(&mut src, &mut buf).unwrap();
        let xml = String::from_utf8(buf).unwrap();

        assert!(xml.contains(r#"<userParam name="openwraw.source_type" value="ESI &amp; APCI"/>"#));
        assert!(xml.contains(r#"<userParam name="acquisition event id" value="42"/>"#));
        assert!(xml.contains(r#"<userParam name="openszraw.segment" value="A &lt; B"/>"#));

        let run_extra = xml
            .find(r#"<userParam name="openwraw.source_type""#)
            .unwrap();
        let spectrum_list = xml.find("<spectrumList ").unwrap();
        assert!(run_extra < spectrum_list);

        let spectrum_start = xml.find("<spectrum ").unwrap();
        let spectrum_xml = &xml[spectrum_start..];
        let last_spectrum_cv = spectrum_xml
            .find(r#"<cvParam cvRef="MS" accession="MS:1000527""#)
            .unwrap();
        let event_id = spectrum_xml
            .find(r#"<userParam name="acquisition event id""#)
            .unwrap();
        let spectrum_extra = spectrum_xml
            .find(r#"<userParam name="openszraw.segment""#)
            .unwrap();
        let scan_list = spectrum_xml.find("<scanList ").unwrap();
        assert!(last_spectrum_cv < event_id);
        assert!(event_id < spectrum_extra);
        assert!(spectrum_extra < scan_list);
    }

    // ---------- Sha1 test vectors -------------------------------------------
    //
    // These vectors come from an independent implementation (Python's
    // `hashlib`), not from this module's own output, so they catch a broken
    // implementation rather than re-confirming it.

    fn sha1_hex(data: &[u8]) -> String {
        let mut h = Sha1::new();
        h.update(data);
        h.finalize().iter().map(|b| format!("{:02x}", b)).collect()
    }

    #[test]
    fn sha1_matches_nist_fips_180_1_vectors() {
        // FIPS 180-1 / RFC 3174 sample vectors.
        assert_eq!(sha1_hex(b"abc"), "a9993e364706816aba3e25717850c26c9cd0d89d");
        assert_eq!(
            sha1_hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
            "84983e441c3bd26ebaae4aa1f95129e5e54670f1"
        );
        let million_a = vec![b'a'; 1_000_000];
        assert_eq!(
            sha1_hex(&million_a),
            "34aa973cd4c4daa4f61eeb2bdbad27316534016f"
        );
    }

    #[test]
    fn sha1_handles_padding_block_boundaries() {
        // Empty input.
        assert_eq!(sha1_hex(b""), "da39a3ee5e6b4b0d3255bfef95601890afd80709");

        // 56 bytes: one byte short of the 64-byte block size, so the 0x80
        // pad byte plus the 8-byte length field must roll over into a
        // second block.
        let msg56: Vec<u8> = (0..56u8).map(|i| b'a' + (i % 26)).collect();
        assert_eq!(sha1_hex(&msg56), "4ad5bb7ae3c4024768d364b77c52128ea3cffebe");

        // 64 bytes: exactly one full block, so padding alone must fill an
        // entire second block.
        let msg64: Vec<u8> = (0..64u8).map(|i| b'a' + (i % 26)).collect();
        assert_eq!(sha1_hex(&msg64), "93249d4c2f8903ebf41ac358473148ae6ddd7042");

        // 1000 bytes: spans several full blocks plus a partial final block.
        let msg1000: Vec<u8> = (0..1000u32).map(|i| (i % 256) as u8).collect();
        assert_eq!(
            sha1_hex(&msg1000),
            "af0b191c2de46fe13fe0908f5a6a4e90e0cafc46"
        );
    }

    #[test]
    fn sha1_streaming_updates_match_single_shot() {
        // CountingWriter feeds Sha1::update() once per `write` call, so the
        // chunk boundaries seen in practice rarely line up with the
        // internal 64-byte block size. Confirm chunking doesn't change the
        // digest.
        let data: Vec<u8> = (0..200u32).map(|i| (i % 256) as u8).collect();

        let mut whole = Sha1::new();
        whole.update(&data);
        let whole_digest = whole.finalize();

        let mut chunked = Sha1::new();
        for chunk in data.chunks(7) {
            chunked.update(chunk);
        }
        let chunked_digest = chunked.finalize();

        assert_eq!(whole_digest, chunked_digest);
    }
}
