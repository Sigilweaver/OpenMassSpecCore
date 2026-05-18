# mass-spec-core

Shared foundation for the open Rust mass-spec parsers
([opentfraw](https://github.com/Sigilweaver/OpenTFRaw),
[opentimstdf](https://github.com/Sigilweaver/OpenTDF),
[openwraw](https://github.com/Sigilweaver/OpenWRaw)).

Defines:

- The vendor-neutral types every parser produces (`SpectrumRecord`,
  `PrecursorInfo`, `ChromatogramRecord`, `RunMetadata`).
- The `SpectrumSource` trait every parser implements.
- One canonical mzML 1.1.0 writer that turns any `SpectrumSource` into a
  valid mzML or indexedmzML document.

Each vendor crate stays a complete standalone tool: a user pulls in
`opentfraw` alone and gets parsing **and** mzML export. This crate is
only the shared vocabulary that keeps the three parsers in lock-step.

See [ROADMAP.md](ROADMAP.md) for the multi-phase plan that this crate is
the foundation of.

## License

Apache-2.0
