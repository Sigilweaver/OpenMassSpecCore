//! Canonical mzML 1.1.0 writer.
//!
//! Consumes any [`SpectrumSource`] and emits valid mzML, with an indexed
//! variant that also writes the `<indexList>` + `<fileChecksum>` trailer so
//! random-access mzML readers (pyteomics, pymzml, sciex-style indexers) work.
//!
//! The implementation here is the vendor-neutral lift of the writer that
//! originally lived in `opentfraw::mzml`. All vendor-specific decisions (the
//! source-file CV term, the native-ID format CV term, the instrument CV term,
//! the per-scan `native_id` strings) come from the source's
//! [`RunMetadata`](crate::RunMetadata) and from each
//! [`SpectrumRecord::native_id`](crate::SpectrumRecord).

use std::io::{Result, Write};

use crate::enums::{Activation, Analyzer, Polarity, ScanMode};
use crate::source::SpectrumSource;
use crate::types::{CvTerm, RunMetadata, SpectrumRecord};

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

fn encode_f64_array(vals: &[f64]) -> String {
    let bytes: Vec<u8> = vals.iter().flat_map(|v| v.to_le_bytes()).collect();
    base64_encode(&bytes)
}

fn encode_f32_array(vals: &[f32]) -> String {
    let bytes: Vec<u8> = vals.iter().flat_map(|v| v.to_le_bytes()).collect();
    base64_encode(&bytes)
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

fn activation_cv(act: Activation, analyzer: Option<Analyzer>) -> (&'static str, &'static str) {
    match act {
        Activation::HCD => ("MS:1000422", "beam-type collision-induced dissociation"),
        Activation::ETD | Activation::EThcD => ("MS:1000598", "electron transfer dissociation"),
        Activation::CID => match analyzer {
            Some(Analyzer::FTMS) => ("MS:1000422", "beam-type collision-induced dissociation"),
            _ => ("MS:1000133", "collision-induced dissociation"),
        },
        Activation::MPID => (
            "MS:1002481",
            "supplemental beam-type collision-induced dissociation",
        ),
        Activation::ECD => ("MS:1000250", "electron capture dissociation"),
        Activation::IRMPD => ("MS:1000262", "infrared multiphoton dissociation"),
        Activation::PD => ("MS:1001880", "in-source collision-induced dissociation"),
        Activation::PQD => ("MS:1000599", "pulsed q dissociation"),
        Activation::UVPD => ("MS:1003246", "ultraviolet photodissociation"),
        Activation::SID => ("MS:1000422", "beam-type collision-induced dissociation"),
    }
}

// ---------- public entry points --------------------------------------------

/// Write the source's spectra as mzML 1.1.0 (un-indexed).
pub fn write_mzml<S: SpectrumSource + ?Sized, W: Write>(src: &mut S, out: &mut W) -> Result<()> {
    let meta = src.run_metadata();
    let count = src.spectrum_count_hint().unwrap_or(0);

    write_prologue(out, &meta, count, false)?;
    for rec in src.iter_spectra() {
        write_spectrum(out, &rec)?;
    }
    writeln!(out, r#"    </spectrumList>"#)?;
    writeln!(out, r#"  </run>"#)?;
    writeln!(out, r#"</mzML>"#)?;
    Ok(())
}

/// Write the source's spectra as indexed mzML 1.1.0 (with `<indexList>` and
/// `<fileChecksum>` trailer).
pub fn write_indexed_mzml<S: SpectrumSource + ?Sized, W: Write>(
    src: &mut S,
    out: &mut W,
) -> Result<()> {
    let meta = src.run_metadata();
    let count = src.spectrum_count_hint().unwrap_or(0);

    let mut cw = CountingWriter::new(out);
    write_prologue(&mut cw, &meta, count, true)?;

    let mut offsets: Vec<(String, u64)> = Vec::with_capacity(count);
    for rec in src.iter_spectra() {
        offsets.push((rec.native_id.clone(), cw.pos));
        write_spectrum(&mut cw, &rec)?;
    }

    writeln!(cw, r#"    </spectrumList>"#)?;
    writeln!(cw, r#"  </run>"#)?;
    writeln!(cw, r#"  </mzML>"#)?;

    let index_list_offset = cw.pos;
    writeln!(cw, r#"  <indexList count="1">"#)?;
    writeln!(cw, r#"    <index name="spectrum">"#)?;
    for (id, offset) in &offsets {
        writeln!(
            cw,
            r#"      <offset idRef="{}">{}</offset>"#,
            escape(id),
            offset
        )?;
    }
    writeln!(cw, r#"    </index>"#)?;
    writeln!(cw, r#"  </indexList>"#)?;

    cw.hashing = false;
    let digest = std::mem::replace(&mut cw.sha1, Sha1::new()).finalize();
    let hex: String = digest.iter().map(|b| format!("{:02x}", b)).collect();

    writeln!(
        cw,
        r#"  <indexListOffset>{}</indexListOffset>"#,
        index_list_offset
    )?;
    writeln!(cw, r#"  <fileChecksum>{}</fileChecksum>"#, hex)?;
    writeln!(cw, r#"</indexedmzML>"#)?;
    Ok(())
}

// ---------- prologue / spectrum body ---------------------------------------

fn write_prologue<W: Write>(
    out: &mut W,
    meta: &RunMetadata,
    n_spectra: usize,
    indexed: bool,
) -> Result<()> {
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
        writeln!(
            out,
            r#"      xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.2_idx.xsd""#
        )?;
        writeln!(out, r#"      version="1.1.0">"#)?;
    }

    writeln!(out, r#"  <cvList count="2">"#)?;
    writeln!(
        out,
        r#"    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="4.1.100" URI="https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"/>"#
    )?;
    writeln!(
        out,
        r#"    <cv id="UO" fullName="Unit Ontology" version="09:04:2014" URI="https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"/>"#
    )?;
    writeln!(out, r#"  </cvList>"#)?;

    writeln!(out, r#"  <fileDescription>"#)?;
    writeln!(out, r#"    <fileContent>"#)?;
    writeln!(
        out,
        r#"      <cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" value=""/>"#
    )?;
    writeln!(
        out,
        r#"      <cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" value=""/>"#
    )?;
    writeln!(out, r#"    </fileContent>"#)?;
    writeln!(out, r#"    <sourceFileList count="1">"#)?;
    writeln!(
        out,
        r#"      <sourceFile id="sf1" name="{}" location="">"#,
        escape(&meta.source_file_name)
    )?;
    write_cv(out, "        ", &meta.source_file_format)?;
    write_cv(out, "        ", &meta.native_id_format)?;
    writeln!(out, r#"      </sourceFile>"#)?;
    writeln!(out, r#"    </sourceFileList>"#)?;
    writeln!(out, r#"  </fileDescription>"#)?;

    writeln!(out, r#"  <softwareList count="1">"#)?;
    writeln!(
        out,
        r#"    <software id="{}" version="{}">"#,
        escape(&meta.software_name),
        escape(&meta.software_version)
    )?;
    writeln!(
        out,
        r#"      <cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="{}"/>"#,
        escape(&meta.software_name)
    )?;
    writeln!(out, r#"    </software>"#)?;
    writeln!(out, r#"  </softwareList>"#)?;

    writeln!(out, r#"  <instrumentConfigurationList count="1">"#)?;
    writeln!(out, r#"    <instrumentConfiguration id="IC1">"#)?;
    write_cv(out, "      ", &meta.instrument)?;
    writeln!(out, r#"    </instrumentConfiguration>"#)?;
    writeln!(out, r#"  </instrumentConfigurationList>"#)?;

    writeln!(out, r#"  <dataProcessingList count="1">"#)?;
    writeln!(out, r#"    <dataProcessing id="dp1">"#)?;
    writeln!(
        out,
        r#"      <processingMethod order="0" softwareRef="{}">"#,
        escape(&meta.software_name)
    )?;
    writeln!(
        out,
        r#"        <cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML" value=""/>"#
    )?;
    writeln!(out, r#"      </processingMethod>"#)?;
    writeln!(out, r#"    </dataProcessing>"#)?;
    writeln!(out, r#"  </dataProcessingList>"#)?;

    writeln!(
        out,
        r#"  <run id="{}" defaultInstrumentConfigurationRef="IC1" defaultSourceFileRef="sf1">"#,
        escape(&meta.source_file_name)
    )?;
    writeln!(
        out,
        r#"    <spectrumList count="{}" defaultDataProcessingRef="dp1">"#,
        n_spectra
    )?;
    Ok(())
}

fn write_cv<W: Write>(out: &mut W, indent: &str, cv: &CvTerm) -> Result<()> {
    writeln!(
        out,
        r#"{indent}<cvParam cvRef="MS" accession="{}" name="{}" value=""/>"#,
        cv.accession,
        escape(&cv.name)
    )
}

fn write_spectrum<W: Write>(out: &mut W, rec: &SpectrumRecord) -> Result<()> {
    let spectrum_type = if rec.ms_level <= 1 {
        ("MS:1000579", "MS1 spectrum")
    } else {
        ("MS:1000580", "MSn spectrum")
    };
    let n_peaks = rec.mz.len();

    writeln!(
        out,
        r#"      <spectrum id="{id}" index="{idx}" defaultArrayLength="{n}">"#,
        id = escape(&rec.native_id),
        idx = rec.index,
        n = n_peaks
    )?;
    writeln!(
        out,
        r#"        <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="{}"/>"#,
        rec.ms_level
    )?;
    writeln!(
        out,
        r#"        <cvParam cvRef="MS" accession="{}" name="{}" value=""/>"#,
        spectrum_type.0, spectrum_type.1
    )?;

    match rec.scan_mode {
        Some(ScanMode::Centroid) => writeln!(
            out,
            r#"        <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>"#
        )?,
        _ => writeln!(
            out,
            r#"        <cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>"#
        )?,
    }

    match rec.polarity {
        Some(Polarity::Positive) => writeln!(
            out,
            r#"        <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" value=""/>"#
        )?,
        Some(Polarity::Negative) => writeln!(
            out,
            r#"        <cvParam cvRef="MS" accession="MS:1000129" name="negative scan" value=""/>"#
        )?,
        None => {}
    }

    let tic = rec.effective_tic();
    let (bp_mz, bp_int) = rec.effective_base_peak().unwrap_or((0.0, 0.0));
    let (lo_mz, hi_mz) = rec.effective_mz_range().unwrap_or((0.0, 0.0));

    writeln!(
        out,
        r#"        <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="{:.6}"/>"#,
        tic
    )?;
    writeln!(
        out,
        r#"        <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="{:.6}"/>"#,
        bp_mz
    )?;
    writeln!(
        out,
        r#"        <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="{:.6}"/>"#,
        bp_int
    )?;
    writeln!(
        out,
        r#"        <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="{:.6}"/>"#,
        lo_mz
    )?;
    writeln!(
        out,
        r#"        <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="{:.6}"/>"#,
        hi_mz
    )?;

    writeln!(out, r#"        <scanList count="1">"#)?;
    writeln!(
        out,
        r#"          <cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>"#
    )?;
    writeln!(out, r#"          <scan>"#)?;

    if let Some(f) = rec.filter.as_deref() {
        if !f.is_empty() {
            writeln!(
                out,
                r#"            <cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="{}"/>"#,
                escape(f)
            )?;
        }
    }

    // mzML stores RT in minutes by convention.
    let rt_min = rec.retention_time_sec / 60.0;
    writeln!(
        out,
        r#"            <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="{:.6}" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>"#,
        rt_min
    )?;

    if let Some(it) = rec.ion_injection_time_ms {
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1000927" name="ion injection time" value="{:.6}" unitCvRef="UO" unitAccession="UO:0000028" unitName="millisecond"/>"#,
            it
        )?;
    }

    if let Some(mob) = rec.inv_mobility {
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1002815" name="inverse reduced ion mobility" value="{:.6}" unitCvRef="MS" unitAccession="MS:1002814" unitName="volt-second per square centimeter"/>"#,
            mob
        )?;
    }

    writeln!(out, r#"            <scanWindowList count="1">"#)?;
    writeln!(out, r#"              <scanWindow>"#)?;
    writeln!(
        out,
        r#"                <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="{:.6}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>"#,
        lo_mz
    )?;
    writeln!(
        out,
        r#"                <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="{:.6}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>"#,
        hi_mz
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

        if pre.target_mz.is_some() || pre.isolation_width.is_some() {
            writeln!(out, r#"            <isolationWindow>"#)?;
            if let Some(mz) = pre.target_mz {
                writeln!(
                    out,
                    r#"              <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{:.6}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>"#,
                    mz
                )?;
            }
            if let Some(w) = pre.isolation_width {
                let half = w / 2.0;
                writeln!(
                    out,
                    r#"              <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="{:.6}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>"#,
                    half
                )?;
                writeln!(
                    out,
                    r#"              <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="{:.6}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>"#,
                    half
                )?;
            }
            writeln!(out, r#"            </isolationWindow>"#)?;
        }

        if let Some(mz) = pre.selected_mz {
            writeln!(out, r#"            <selectedIonList count="1">"#)?;
            writeln!(out, r#"              <selectedIon>"#)?;
            writeln!(
                out,
                r#"                <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="{:.6}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>"#,
                mz
            )?;
            if let Some(z) = pre.charge {
                writeln!(
                    out,
                    r#"                <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="{z}"/>"#
                )?;
            }
            if let Some(i) = pre.intensity {
                writeln!(
                    out,
                    r#"                <cvParam cvRef="MS" accession="MS:1000042" name="peak intensity" value="{:.6}"/>"#,
                    i
                )?;
            }
            writeln!(out, r#"              </selectedIon>"#)?;
            writeln!(out, r#"            </selectedIonList>"#)?;
        }

        writeln!(out, r#"            <activation>"#)?;
        if let Some(act) = pre.activation {
            let (acc, name) = activation_cv(act, pre.analyzer);
            writeln!(
                out,
                r#"              <cvParam cvRef="MS" accession="{acc}" name="{name}" value=""/>"#
            )?;
        } else {
            writeln!(
                out,
                r#"              <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" value=""/>"#
            )?;
        }
        if let Some(e) = pre.collision_energy {
            if pre.ce_is_nce {
                writeln!(
                    out,
                    r#"              <cvParam cvRef="MS" accession="MS:1002013" name="normalized collision energy" value="{:.2}"/>"#,
                    e
                )?;
            } else {
                writeln!(
                    out,
                    r#"              <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="{:.2}" unitCvRef="UO" unitAccession="UO:0000266" unitName="electronvolt"/>"#,
                    e
                )?;
            }
        }
        writeln!(out, r#"            </activation>"#)?;
        writeln!(out, r#"          </precursor>"#)?;
        writeln!(out, r#"        </precursorList>"#)?;
    }

    if n_peaks > 0 {
        let mz_b64 = encode_f64_array(&rec.mz);
        let int_b64 = encode_f32_array(&rec.intensity);

        writeln!(out, r#"        <binaryDataArrayList count="2">"#)?;

        writeln!(
            out,
            r#"          <binaryDataArray encodedLength="{}">"#,
            mz_b64.len()
        )?;
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value=""/>"#
        )?;
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>"#
        )?;
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>"#
        )?;
        writeln!(out, r#"            <binary>{mz_b64}</binary>"#)?;
        writeln!(out, r#"          </binaryDataArray>"#)?;

        writeln!(
            out,
            r#"          <binaryDataArray encodedLength="{}">"#,
            int_b64.len()
        )?;
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value=""/>"#
        )?;
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>"#
        )?;
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>"#
        )?;
        writeln!(out, r#"            <binary>{int_b64}</binary>"#)?;
        writeln!(out, r#"          </binaryDataArray>"#)?;

        writeln!(out, r#"        </binaryDataArrayList>"#)?;
    }

    writeln!(out, r#"      </spectrum>"#)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base64_rfc_vectors() {
        assert_eq!(base64_encode(b""), "");
        assert_eq!(base64_encode(b"f"), "Zg==");
        assert_eq!(base64_encode(b"fo"), "Zm8=");
        assert_eq!(base64_encode(b"foo"), "Zm9v");
        assert_eq!(base64_encode(b"foobar"), "Zm9vYmFy");
        assert_eq!(base64_encode(b"Man"), "TWFu");
    }
}
