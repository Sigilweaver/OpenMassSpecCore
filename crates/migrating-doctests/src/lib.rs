//! Compiles the `rust` blocks of `MIGRATING-4.0.md` as doc tests, so the
//! "After" snippets stay true to the API. This crate is not published; it
//! lives beside the published crates so that no published crate references a
//! file outside its own package.

#![doc = include_str!("../../../MIGRATING-4.0.md")]
