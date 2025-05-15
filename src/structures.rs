// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! This file contains public data structures used in the `minitpr` library.

pub use mendeleev::Element;

use crate::DIM;

/// Structure representing the TPR file.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TprFile {
    /// TPR file header.
    pub header: TprHeader,
    /// Name of the molecular system.
    pub system_name: String,
    /// Dimensions of the simulation box.
    pub simbox: Option<SimBox>,
    /// System topology.
    pub topology: TprTopology,
}

/// Structure representing the header of the TPR file.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TprHeader {
    /// Gromacs version used to write the tpr file.
    pub gromacs_version: String,
    /// Precision of the tpr file.
    pub precision: Precision,
    /// Version of the tpr file.
    pub tpr_version: i32,
    /// Generation of the tpr file.
    pub tpr_generation: i32,
    /// Tpr file tag.
    pub file_tag: String,
    /// Number of atoms.
    pub n_atoms: i32,
    /// Number of temperature coupling groups.
    pub n_coupling_groups: i32,
    /// Value of the alchemical state.
    pub fep_state: i32,
    /// Value of lambda.
    pub lambda: f64,
    /// Is input record present?
    pub has_input_record: bool,
    /// Is topology present?
    pub has_topology: bool,
    /// Are positions present?
    pub has_positions: bool,
    /// Are velocities present?
    pub has_velocities: bool,
    /// Are forces present?
    pub has_forces: bool,
    /// Is the simulation box present?
    pub has_box: bool,
    /// Size of the body of the tpr file.
    pub body_size: Option<i64>,
}

/// Structure representing the topology of the TPR file.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TprTopology {
    /// List of atoms in the system.
    pub atoms: Vec<Atom>,
    /// List of bonds between atoms in the system.
    /// The order of bonds is undefined.
    pub bonds: Vec<Bond>,
}

/// Structure representing simulation box dimensions.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SimBox {
    pub simbox: [[f64; DIM]; DIM],
    pub simbox_rel: [[f64; DIM]; DIM],
    pub simbox_v: [[f64; DIM]; DIM],
}

/// Enum representing precision of the tpr file.
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Precision {
    Single,
    Double,
}

/// Structure representing an atom.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Atom {
    /// Name of the atom.
    pub atom_name: String,
    /// Atom number. All atoms are numbered sequentially, starting from 1.
    pub atom_number: i32,
    /// Name of the residue this atom is part of.
    pub residue_name: String,
    /// Residue number. All residues are numbered sequentially, starting from 1.
    pub residue_number: i32,
    /// Mass of the atom.
    pub mass: f64,
    /// Charge of the atom.
    pub charge: f64,
    /// Element this atom belongs to.
    pub element: Option<Element>,
    /// Position of the atom.
    pub position: Option<[f64; 3]>,
    /// Velocity of the atom.
    pub velocity: Option<[f64; 3]>,
    /// Force acting on the atom.
    pub force: Option<[f64; 3]>,
}

/// Structure representing a bond between atoms.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Bond {
    /// Global index of the first atom involved in the bond.
    /// Note the word is `index`, not `number`! Atom indices start from 0.
    /// Indices correspond to indices of the atoms in the `TprTopology::atoms` vector.
    pub atom1: usize,
    /// Global index of the second atom involved in the bond.
    pub atom2: usize,
}
