//! Integration coverage for the indexed-mzML `<fileChecksum>` trailer
//! (issue #11).
//!
//! `write_indexed_mzml` computes its `fileChecksum` with a hand-rolled,
//! internal SHA-1 (see `src/mzml.rs`) so the crate doesn't need a crypto
//! dependency just for this one footer. That implementation has its own
//! unit tests against NIST vectors in `src/mzml.rs`, but nothing exercised
//! the actual code path that feeds a generated document's bytes through it.
//! This test closes that gap: it generates a small indexed-mzML document
//! through the public API, recomputes SHA-1 over the same bytes with a
//! second, independent implementation written here, and asserts the two
//! agree.
//!
//! Note: the current implementation hashes everything written up to and
//! including the closing `</indexList>` tag, and stops there - it does
//! not include the `<indexListOffset>` element or the `<fileChecksum>`
//! open tag in the digest. This test locks down that actual behavior (so
//! a refactor can't silently change it without a test failing) rather
//! than asserting what the bundled schema doc-comment describes as the
//! nominal scope; see the accompanying PR discussion for that discrepancy.

use openmassspec_core::{
    write_indexed_mzml, CvTerm, Polarity, RunMetadata, ScanMode, SpectrumRecord, SpectrumSource,
};

// ---------- independent SHA-1 (whole-message, non-streaming) ---------------
//
// Deliberately not shared with `src/mzml.rs`'s streaming implementation:
// the point of this test is to catch a bug in *that* implementation, which
// it can't do if it reuses the same code to check itself.

fn sha1_hex(message: &[u8]) -> String {
    let mut h0: u32 = 0x67452301;
    let mut h1: u32 = 0xEFCDAB89;
    let mut h2: u32 = 0x98BADCFE;
    let mut h3: u32 = 0x10325476;
    let mut h4: u32 = 0xC3D2E1F0;

    let bit_len = (message.len() as u64) * 8;
    let mut padded = message.to_vec();
    padded.push(0x80);
    while padded.len() % 64 != 56 {
        padded.push(0);
    }
    padded.extend_from_slice(&bit_len.to_be_bytes());

    for block in padded.chunks(64) {
        let mut w = [0u32; 80];
        for (i, word) in w.iter_mut().enumerate().take(16) {
            *word = u32::from_be_bytes(block[i * 4..i * 4 + 4].try_into().unwrap());
        }
        for i in 16..80 {
            w[i] = (w[i - 3] ^ w[i - 8] ^ w[i - 14] ^ w[i - 16]).rotate_left(1);
        }

        let (mut a, mut b, mut c, mut d, mut e) = (h0, h1, h2, h3, h4);
        for (i, &wi) in w.iter().enumerate() {
            let (f, k) = match i {
                0..=19 => ((b & c) | ((!b) & d), 0x5A827999u32),
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

        h0 = h0.wrapping_add(a);
        h1 = h1.wrapping_add(b);
        h2 = h2.wrapping_add(c);
        h3 = h3.wrapping_add(d);
        h4 = h4.wrapping_add(e);
    }

    [h0, h1, h2, h3, h4]
        .iter()
        .map(|w| format!("{:08x}", w))
        .collect()
}

// ---------- minimal SpectrumSource for a one-spectrum document -------------

struct OneSpectrumSource {
    meta: RunMetadata,
    spectrum: SpectrumRecord,
}

impl OneSpectrumSource {
    fn new() -> Self {
        Self {
            meta: RunMetadata {
                source_file_name: "checksum-fixture.raw".into(),
                source_file_format: CvTerm::new("MS:1000563", "Thermo RAW format"),
                native_id_format: CvTerm::new("MS:1000768", "Thermo nativeID format"),
                instrument: CvTerm::new("MS:1001911", "Q Exactive"),
                instrument_serial_number: None,
                acquisition_software_name: None,
                acquisition_software_version: None,
                software_name: "checksum-test".into(),
                software_version: "0.0.0".into(),
                start_timestamp: None,
                mobility_array_kind: None,
                analyzers: Vec::new(),
            },
            spectrum: SpectrumRecord {
                index: 0,
                scan_number: 1,
                native_id: "controllerType=0 controllerNumber=1 scan=1".into(),
                ms_level: 1,
                polarity: Some(Polarity::Positive),
                scan_mode: Some(ScanMode::Centroid),
                analyzer: None,
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
                mz: vec![100.0, 200.0, 300.0],
                intensity: vec![10.0, 5.0, 1.0],
                inv_mobility_per_peak: None,
            },
        }
    }
}

impl SpectrumSource for OneSpectrumSource {
    fn run_metadata(&self) -> RunMetadata {
        self.meta.clone()
    }

    fn iter_spectra<'a>(&'a mut self) -> Box<dyn Iterator<Item = SpectrumRecord> + 'a> {
        Box::new(std::iter::once(self.spectrum.clone()))
    }

    fn spectrum_count_hint(&self) -> Option<usize> {
        Some(1)
    }
}

#[test]
fn embedded_file_checksum_matches_independent_sha1_recompute() {
    let mut src = OneSpectrumSource::new();
    let mut buf = Vec::new();
    write_indexed_mzml(&mut src, &mut buf).unwrap();

    let needle = b"</indexList>\n";
    let end_of_index_list = buf
        .windows(needle.len())
        .position(|w| w == needle)
        .expect("generated document must contain a closing </indexList> tag")
        + needle.len();

    let expected = sha1_hex(&buf[..end_of_index_list]);

    let s = String::from_utf8(buf.clone()).unwrap();
    let open = s
        .find("<fileChecksum>")
        .expect("fileChecksum element present")
        + "<fileChecksum>".len();
    let close = s[open..]
        .find("</fileChecksum>")
        .expect("fileChecksum closing tag present")
        + open;
    let embedded = &s[open..close];

    assert_eq!(
        embedded, expected,
        "embedded <fileChecksum> should equal an independently-recomputed \
         SHA-1 of the bytes written up to and including </indexList>"
    );
    // Sanity: a real 40-char lowercase hex SHA-1 digest, not an empty or
    // truncated string.
    assert_eq!(embedded.len(), 40);
    assert!(embedded.chars().all(|c| c.is_ascii_hexdigit()));
}
