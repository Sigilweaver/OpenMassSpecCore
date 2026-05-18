//! Aggregate error type for the OpenProteo stack.
//!
//! Each vendor crate defines its own narrow `Error` enum. When code needs to
//! handle errors from multiple vendors uniformly - the umbrella `openproteo-io`
//! crate, the `vendor2mzml` CLI, ProLance ingest - it converts those into
//! [`Error`] (this aggregate). Downstream users get a single error vocabulary
//! and `?`-propagates cleanly across vendor boundaries.

use std::io;

/// Stack-wide error type. Vendor errors are erased into a boxed trait object
/// so this enum stays version-stable as vendor crates evolve.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// I/O failure (file not found, permission denied, short read, ...).
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),

    /// Vendor parser failed. The inner error is the vendor crate's own error
    /// type, preserved via `Box<dyn std::error::Error>` so callers can
    /// downcast if they need vendor-specific context.
    #[error("vendor parser error: {0}")]
    Vendor(Box<dyn std::error::Error + Send + Sync + 'static>),

    /// Format detection failed or the input does not match a supported vendor.
    #[error("unsupported or unrecognized format: {0}")]
    Format(String),

    /// Conformance harness rejected a record stream.
    #[error("conformance violation: {0}")]
    Conformance(String),
}

/// Convenience alias.
pub type Result<T> = std::result::Result<T, Error>;

impl Error {
    /// Wrap a vendor-crate error into [`Error::Vendor`].
    pub fn vendor<E>(err: E) -> Self
    where
        E: std::error::Error + Send + Sync + 'static,
    {
        Error::Vendor(Box::new(err))
    }
}
