// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

//! This file contains functions for obtaining molecule blocks from TPR file.

use crate::{
    errors::ParseTprError,
    structures::{Atom, Bond, Precision},
};

use super::{moltypes::MoleculeType, xdr::XdrFile};

/// Structure representing a molecule block.
#[derive(Debug, Clone)]
pub(super) struct MolBlock {
    pub molecule_type: i32,
    pub n_molecules: i32,
}

impl MolBlock {
    /// Get `MolBlock` from `XdrFile`.
    pub(super) fn parse(
        xdrfile: &mut XdrFile,
        precision: Precision,
    ) -> Result<Self, ParseTprError> {
        let molecule_type = xdrfile.read_i32()?;
        let n_molecules = xdrfile.read_i32()?;
        // ignore n_atoms_per_molecule field
        xdrfile.jump(4)?;

        // skip position restraints
        for _ in 0..2 {
            let n_posres = xdrfile.read_i32()?;
            xdrfile.skip_multiple_reals(precision, crate::DIM as i64 * n_posres as i64)?;
        }

        Ok(MolBlock {
            molecule_type,
            n_molecules,
        })
    }

    /// Unpack `MolBlock` to molecules, i.e., a vector of atoms and a vector of bonds.
    pub(super) fn unpack2molecules(
        &self,
        molecule_types: &[MoleculeType],
        atom_counter: &mut i32,
        residue_counter: &mut i32,
    ) -> Result<(Vec<Atom>, Vec<Bond>), ParseTprError> {
        let moltype = match molecule_types.get(self.molecule_type as usize) {
            Some(x) => x,
            None => return Err(ParseTprError::CouldNotConstructTopology),
        };

        let mut atoms = Vec::with_capacity(moltype.atoms.len() * self.n_molecules as usize);
        let mut bonds = Vec::new();

        for _ in 0..self.n_molecules {
            let (new_atoms, new_bonds) = moltype.unpack2molecule(atom_counter, residue_counter)?;
            atoms.extend(new_atoms);
            bonds.extend(new_bonds);
        }

        Ok((atoms, bonds))
    }
}
