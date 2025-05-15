// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! # minitpr
//!
//! `minitpr` is a pure Rust library designed for parsing Gromacs tpr files, focusing on extracting system topology and structure.
//!
//! ## Capabilities and Limitations
//! - Supports parsing of tpr files from version 103 onwards (Gromacs 5.1 and later).
//! - Extracts system topology and structure: atoms, their basic properties (including positions, velocities, and forces), and bonds between atoms (including intermolecular bonds).
//! - Does **not** support parsing of force-field and simulation parameters, nor does it offer capabilities to write tpr files.
//!
//! ## Usage
//!
//! To include `minitpr` in your project, add it as a dependency using Cargo:
//!
//! ```shell
//! cargo add minitpr
//!```
//!
//! ### Parsing a tpr file
//! Example usage to parse a tpr file and handle potential errors:
//! ```rust
//! use minitpr::TprFile;
//!
//! fn main() {
//!     let tpr = match TprFile::parse("topol.tpr") {
//!         Ok(file) => file,
//!         Err(error) => {
//!             eprintln!("{}", error);
//!             return;
//!         },
//!     };
//!
//!     // now you can work with the `tpr` object, accessing its properties and data
//!     // for instance, to iterate through the atoms of the molecular system, use:
//!     for atom in tpr.topology.atoms.iter() {
//!         // perform some operation with the atom
//!     }
//! }
//! ```
//!
//! ### Data Structures
//! The [`TprFile`](`crate::TprFile`) structure encapsulates the following information:
//!
//! - Header: Metadata about the tpr file (see [`TprHeader`](`crate::TprHeader`) structure).
//! - Molecular System Name: The name of the simulated system.
//! - Simulation Box Dimensions: Available within the [`SimBox`](`crate::SimBox`) structure if present.
//! - System Topology: Topology of the molecular system containing atoms and bonds (see [`TprTopology`](`crate::TprTopology`) structure).
//!
//! Each atom (see [`Atom`](`crate::Atom`)) represented in the system topology includes:
//! - Atom name.
//! - Sequential atom number, starting from 1.
//! - Residue name.
//! - Sequential residue number, starting from 1.
//! - Mass.
//! - Charge.
//! - Element (`None` if unidentifiable).
//! - Position (`None` if not present).
//! - Velocity (`None` if not present).
//! - Force (`None` if not present).
//!
//! ## Features
//! ### Serialization/Deserialization
//! Enable (de)serialization support for `TprFile` with `serde` by adding the feature flag during installation:
//! ```shell
//! cargo add minitpr --features serde
//! ```
//!
//! ## License
//! `minitpr` is open-sourced under either the [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0) or the [MIT License](https://opensource.org/license/MIT) at your option.
//!
//! ## Disclaimer
//! Due to the huge number of various force field and simulation parameters and Gromacs options, it is very difficult to *comprehensively* test the `minitpr` library.
//! Your contributions in the form of tests, tpr files with unusual parameters, or new functionality are very welcome. See [github.com/Ladme/minitpr](https://github.com/Ladme/minitpr).
//!
//! If the library is unable to parse your tpr file, but you believe it should be able to, please open a [GitHub issue](https://github.com/Ladme/minitpr/issues) and **upload your tpr file**.
//!

use errors::ParseTprError;
use std::path::Path;

pub mod errors;
mod parse;
pub mod structures;

pub use structures::*;

/// Current version of the `minitpr` library.
pub const MINITPR_VERSION: &str = env!("CARGO_PKG_VERSION");

/// Number of spatial dimensions.
pub(crate) const DIM: usize = 3;
/// Number of fields in the `F_RBDIHS` and `F_FOURDIHS` function types
pub(crate) const NR_RBDIHS: usize = 6;
/// Number of fields in the `F_CBTDIHS` function type
pub(crate) const NR_CBTDIHS: usize = 6;
/// Number of group types (TemperatureCoupling, EnergyOutput, Acceleration, etc.).
pub(crate) const NR_GROUP_TYPES: usize = 10;

impl TprFile {
    /// Parse a Gromacs tpr file.
    ///
    /// ## Parameters
    /// - `filename`: path to the tpr file to read
    ///
    /// ## Returns
    /// - [`TprFile`](`crate::TprFile`) structure, if successful.
    /// - Otherwise [`ParseTprError`](`crate::errors::ParseTprError`).
    ///
    /// ## Notes
    /// - Only tpr files version 103 or higher are supported (Gromacs 5.1 onwards).
    /// - The function only parses the following information: tpr file header,
    ///   name of the system, simulation box, system topology (atoms and bonds),
    ///   positions, velocities, and forces (if present).
    /// - Force-field properties and simulation parameters are NOT parsed.
    /// - If the tpr file does not contain topology information, this function will return an error.
    pub fn parse(filename: impl AsRef<Path>) -> Result<Self, ParseTprError> {
        parse::parse_tpr(filename)
    }
}
