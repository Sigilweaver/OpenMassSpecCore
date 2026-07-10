# Security Policy

## Supported Versions

| Version | Supported |
| ------- | --------- |
| latest  | Yes       |
| older   | No        |

Only the latest published release receives security updates.

## Reporting a Vulnerability

**Do not open a public GitHub issue for security vulnerabilities.**

Report privately via [GitHub Security Advisories](https://github.com/Sigilweaver/OpenMassSpecCore/security/advisories/new).

Include:

- A description of the vulnerability and its potential impact.
- Steps to reproduce or a proof of concept.
- The OS, Rust toolchain, and crate version you were running.

Expect an initial acknowledgment within 7 days.

## Scope

In scope:

- **mzML writer correctness on adversarial input.** `openmassspec-core`
  ships the canonical mzML 1.1.0 writer that every vendor parser
  funnels through. Crashes (panics, OOB writes), arbitrary file
  writes, or XML-injection issues triggered by crafted
  `SpectrumRecord` data are in scope.
- **Arrow bridge soundness** (under the `arrow` feature): unsound
  memory access, undefined behavior, or buffer-length confusion when
  building `RecordBatch` values from spectra.
- **Conformance harness correctness**: false-positive *passes* (the
  harness lets invalid input through). False negatives are bugs but
  not security issues.
- **Supply-chain integrity** of published artifacts on crates.io.

Out of scope:

- Vulnerabilities in third-party crates with no demonstrated exploit
  path through this crate. Forward those upstream.
- Denial of service via legitimately oversized spectrum streams.
- Issues that require write access to the source tree.

## Disclosure

We follow coordinated disclosure. Reporters are credited in the
release notes unless they prefer to remain anonymous. We aim to ship
a fix within 30 days of confirming a high or critical issue.

## Stack context

`openmassspec-core` is the foundation crate of the
[OpenMassSpec](https://github.com/Sigilweaver/OpenMassSpec) stack. Reports
that involve vendor-specific parsing are usually better routed to the
relevant parser repo:

- Thermo `.raw`: [OpenTFRaw](https://github.com/Sigilweaver/OpenTFRaw)
- Bruker `.d/` (timsTOF): [OpenTimsTDF](https://github.com/Sigilweaver/OpenTimsTDF)
- Waters MassLynx `.raw/`: [OpenWRaw](https://github.com/Sigilweaver/OpenWRaw)
