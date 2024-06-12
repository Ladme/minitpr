// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file defines errors that can be returned by the `minitpr` library.

use std::path::Path;

use colored::{ColoredString, Colorize};
use thiserror::Error;

fn path_to_yellow(path: &Path) -> ColoredString {
    path.to_str().unwrap().yellow()
}

/// Errors that can occur when parsing a tpr file.
#[derive(Error, Debug)]
pub enum ParseTprError {
    /// Used when the tpr file could not be opened.
    #[error("{} file '{}' could not be opened for reading", "error:".red().bold(), path_to_yellow(.0))]
    CouldNotOpen(Box<Path>),
    /// Used when expected data could not be read from a tpr file.
    #[error("{} could not read data from a tpr file (`{}`)", "error:".red().bold(), .0.to_string().yellow())]
    CouldNotRead(#[from] std::io::Error),
    /// Used when the file is not a tpr file.
    #[error("{} parsed file is not a tpr file", "error:".red().bold())]
    NotTpr,
    /// Used when the precision of the tpr file is not supported.
    #[error("{} unsupported tpr file precision `{}`", "error:".red().bold(), .0.to_string().yellow())]
    UnsupportedPrecision(i32),
    /// Used when the version of the tpr file is not supported (is older than version 103).
    #[error("{} unsupported tpr file version `{}`", "error:".red().bold(), .0.to_string().yellow())]
    UnsupportedVersion(i32),
    /// Used when a symbol is requested from the SymTable that does not exist.
    #[error("{} invalid SymTable call: `{}` is out-of-range of the SymTable", "error:".red().bold(), .0.to_string().yellow())]
    IndexNotInSymTable(i32),
    /// Used when sanity check for Interaction parsing fails.
    #[error("{} discrepancy in Interaction of type `{}`: the number of instances is not divisible by the number of interacting atoms + 1",
    "error:".red().bold(), .0.to_string().yellow())]
    InteractionDiscrepancy(i32),
    /// Used when `interaction_type_index` for a Interaction does not exist.
    #[error("{} interaction type index `{}` does not exist", "error:".red().bold(), .0.to_string().yellow())]
    InvalidInteractionType(i32),
    /// Used when the tpr file has been parsed seemingly successfully but topology could not be constructed.
    #[error("{} could not construct molecular topology", "error:".red().bold())]
    CouldNotConstructTopology,
    /// Used when there is an inconsistency in the number of atoms read from the TPR file.
    #[error("{} inconsistent number of atoms in the tpr file (expected `{}` atoms, got `{}` atoms)", "error:".red().bold(), .0.to_string().yellow(), .1.to_string().yellow())]
    InconsistentNumberOfAtoms(i32, i32),
    /// Used when an interaction classified as `bond` is involving different number of atoms than 2.
    #[error("{} invalid number of atoms (`{}`) involved in a bond", "error:".red().bold(), .0.to_string().yellow())]
    InvalidNumberOfBondedAtoms(usize),
    /// Used when the size of intermolecular exclusion group is negative.
    #[error("{} invalid intermolecular exclusion group size (expected a positive value, got `{}`)", "error:".red().bold(), .0.to_string().yellow())]
    InvalidIntermolecularExclusionGroupSize(i64),
}
