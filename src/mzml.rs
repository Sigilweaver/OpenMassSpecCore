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

use crate::enums::{Activation, Analyzer, MobilityArrayKind, Polarity, ScanMode};
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
    let mobility_kind = meta.mobility_array_kind;

    write_prologue(out, &meta, count, false)?;
    for rec in src.iter_spectra() {
        write_spectrum(out, &rec, mobility_kind)?;
    }
    writeln!(out, r#"    </spectrumList>"#)?;

    let chroms: Vec<ChromatogramRecord> = src.iter_chromatograms().collect();
    if !chroms.is_empty() {
        writeln!(
            out,
            r#"    <chromatogramList count="{}" defaultDataProcessingRef="dp1">"#,
            chroms.len()
        )?;
        for rec in &chroms {
            write_chromatogram(out, rec)?;
        }
        writeln!(out, r#"    </chromatogramList>"#)?;
    }

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
    let mobility_kind = meta.mobility_array_kind;

    let mut cw = CountingWriter::new(out);
    write_prologue(&mut cw, &meta, count, true)?;

    let mut offsets: Vec<(String, u64)> = Vec::with_capacity(count);
    for rec in src.iter_spectra() {
        offsets.push((rec.native_id.clone(), cw.pos));
        write_spectrum(&mut cw, &rec, mobility_kind)?;
    }

    writeln!(cw, r#"    </spectrumList>"#)?;

    let chroms: Vec<ChromatogramRecord> = src.iter_chromatograms().collect();
    let mut chrom_offsets: Vec<(String, u64)> = Vec::with_capacity(chroms.len());
    if !chroms.is_empty() {
        writeln!(
            cw,
            r#"    <chromatogramList count="{}" defaultDataProcessingRef="dp1">"#,
            chroms.len()
        )?;
        for rec in &chroms {
            chrom_offsets.push((rec.id.clone(), cw.pos));
            write_chromatogram(&mut cw, rec)?;
        }
        writeln!(cw, r#"    </chromatogramList>"#)?;
    }

    writeln!(cw, r#"  </run>"#)?;
    writeln!(cw, r#"  </mzML>"#)?;

    let index_list_offset = cw.pos;
    let n_indexes = 1 + usize::from(!chrom_offsets.is_empty());
    writeln!(cw, r#"  <indexList count="{n_indexes}">"#)?;
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
    if !chrom_offsets.is_empty() {
        writeln!(cw, r#"    <index name="chromatogram">"#)?;
        for (id, offset) in &chrom_offsets {
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
        // Non-indexed documents have a plain <mzML> root, so they reference
        // the base schema (matching version="1.1.0"), not the indexed
        // wrapper schema used for <indexedmzML> output.
        writeln!(
            out,
            r#"      xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd""#
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

    match &meta.start_timestamp {
        Some(ts) => writeln!(
            out,
            r#"  <run id="{}" defaultInstrumentConfigurationRef="IC1" defaultSourceFileRef="sf1" startTimeStamp="{}">"#,
            escape(&meta.source_file_name),
            escape(ts)
        )?,
        None => writeln!(
            out,
            r#"  <run id="{}" defaultInstrumentConfigurationRef="IC1" defaultSourceFileRef="sf1">"#,
            escape(&meta.source_file_name)
        )?,
    }
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

fn write_spectrum<W: Write>(
    out: &mut W,
    rec: &SpectrumRecord,
    mobility_kind: Option<MobilityArrayKind>,
) -> Result<()> {
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

    if let Some(cv) = rec.faims_cv {
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1001581" name="FAIMS compensation voltage" value="{:.4}" unitCvRef="UO" unitAccession="UO:0000218" unitName="volt"/>"#,
            cv
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
        let mobility_b64_opt = rec
            .inv_mobility_per_peak
            .as_ref()
            .filter(|v| v.len() == n_peaks)
            .map(|v| encode_f32_array(v));
        let array_count = 2 + usize::from(mobility_b64_opt.is_some());

        writeln!(
            out,
            r#"        <binaryDataArrayList count="{array_count}">"#
        )?;

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

        if let Some(mobility_b64) = mobility_b64_opt {
            let (cv_acc, cv_name, unit_acc, unit_ref, unit_name) = match mobility_kind {
                Some(MobilityArrayKind::DriftTimeMilliseconds) => (
                    "MS:1003007",
                    "raw ion mobility array",
                    "UO:0000028",
                    "UO",
                    "millisecond",
                ),
                Some(MobilityArrayKind::InverseReducedVsPerCm2) | None => (
                    "MS:1003008",
                    "raw inverse reduced ion mobility array",
                    "MS:1002814",
                    "MS",
                    "volt-second per square centimeter",
                ),
            };
            writeln!(
                out,
                r#"          <binaryDataArray encodedLength="{}">"#,
                mobility_b64.len()
            )?;
            writeln!(
                out,
                r#"            <cvParam cvRef="MS" accession="{cv_acc}" name="{cv_name}" value="" unitCvRef="{unit_ref}" unitAccession="{unit_acc}" unitName="{unit_name}"/>"#
            )?;
            writeln!(
                out,
                r#"            <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>"#
            )?;
            writeln!(
                out,
                r#"            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>"#
            )?;
            writeln!(out, r#"            <binary>{mobility_b64}</binary>"#)?;
            writeln!(out, r#"          </binaryDataArray>"#)?;
        }

        writeln!(out, r#"        </binaryDataArrayList>"#)?;
    }

    writeln!(out, r#"      </spectrum>"#)?;
    Ok(())
}

fn write_chromatogram<W: Write>(out: &mut W, rec: &ChromatogramRecord) -> Result<()> {
    let n = rec.time_sec.len();

    writeln!(
        out,
        r#"      <chromatogram id="{id}" index="{idx}" defaultArrayLength="{n}">"#,
        id = escape(&rec.id),
        idx = rec.index,
    )?;

    match &rec.chromatogram_type {
        Some(cv) => write_cv(out, "        ", cv)?,
        None => writeln!(
            out,
            r#"        <cvParam cvRef="MS" accession="MS:1000235" name="total ion current chromatogram" value=""/>"#
        )?,
    }

    if let Some(mz) = rec.precursor_mz {
        writeln!(out, r#"        <precursor>"#)?;
        writeln!(out, r#"          <isolationWindow>"#)?;
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{:.6}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>"#,
            mz
        )?;
        writeln!(out, r#"          </isolationWindow>"#)?;
        writeln!(out, r#"          <activation>"#)?;
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" value=""/>"#
        )?;
        writeln!(out, r#"          </activation>"#)?;
        writeln!(out, r#"        </precursor>"#)?;
    }

    if let Some(mz) = rec.product_mz {
        writeln!(out, r#"        <product>"#)?;
        writeln!(out, r#"          <isolationWindow>"#)?;
        writeln!(
            out,
            r#"            <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{:.6}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>"#,
            mz
        )?;
        writeln!(out, r#"          </isolationWindow>"#)?;
        writeln!(out, r#"        </product>"#)?;
    }

    // binaryDataArrayList is required (not optional) on ChromatogramType,
    // unlike SpectrumType where it's dropped for zero-peak spectra.
    // mzML stores time arrays in minutes by convention, matching scan start
    // time; `time_sec` is stored in seconds per `ChromatogramRecord`'s docs.
    let time_min: Vec<f32> = rec.time_sec.iter().map(|&t| t / 60.0).collect();
    let time_b64 = encode_f32_array(&time_min);
    let int_b64 = encode_f32_array(&rec.intensity);

    writeln!(out, r#"        <binaryDataArrayList count="2">"#)?;

    writeln!(
        out,
        r#"          <binaryDataArray encodedLength="{}">"#,
        time_b64.len()
    )?;
    writeln!(
        out,
        r#"            <cvParam cvRef="MS" accession="MS:1000595" name="time array" value="" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>"#
    )?;
    writeln!(
        out,
        r#"            <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" value=""/>"#
    )?;
    writeln!(
        out,
        r#"            <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>"#
    )?;
    writeln!(out, r#"            <binary>{time_b64}</binary>"#)?;
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
        start_timestamp: Option<String>,
        spectra: Vec<SpectrumRecord>,
        chroms: Vec<ChromatogramRecord>,
    }

    impl ToySource {
        fn new() -> Self {
            Self {
                start_timestamp: None,
                spectra: vec![minimal_spectrum(0, None)],
                chroms: Vec::new(),
            }
        }
    }

    impl SpectrumSource for ToySource {
        fn run_metadata(&self) -> RunMetadata {
            RunMetadata {
                source_file_name: "toy.raw".into(),
                source_file_format: CvTerm::new("MS:1000563", "Thermo RAW format"),
                native_id_format: CvTerm::new("MS:1000768", "Thermo nativeID format"),
                instrument: CvTerm::new("MS:1001911", "Q Exactive"),
                software_name: "toy".into(),
                software_version: "0.0.0".into(),
                start_timestamp: self.start_timestamp.clone(),
                mobility_array_kind: None,
            }
        }

        fn iter_spectra<'a>(&'a mut self) -> Box<dyn Iterator<Item = SpectrumRecord> + 'a> {
            Box::new(self.spectra.clone().into_iter())
        }

        fn iter_chromatograms<'a>(
            &'a mut self,
        ) -> Box<dyn Iterator<Item = ChromatogramRecord> + 'a> {
            Box::new(self.chroms.clone().into_iter())
        }

        fn spectrum_count_hint(&self) -> Option<usize> {
            Some(self.spectra.len())
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
}
