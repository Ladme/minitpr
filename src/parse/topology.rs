// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file contains functions for obtaining system topology from a TPR file.

use super::{ffparams::FFParams, molblocks::MolBlock, moltypes::MoleculeType, xdr::XdrFile};
use crate::{
    errors::ParseTprError,
    structures::{Precision, TprTopology},
};

use super::symtab::SymTable;

impl TprTopology {
    /// Get system topology from the tpr file.
    pub(super) fn parse(
        xdrfile: &mut XdrFile,
        precision: Precision,
        tpr_version: i32,
        symbol_table: &SymTable,
        ffparams: &FFParams,
        expected_n_atoms: i32,
    ) -> Result<TprTopology, ParseTprError> {
        // get molecule types
        let n_moltypes = xdrfile.read_i32()?;

        let mut molecule_types = Vec::with_capacity(n_moltypes as usize);
        for _ in 0..n_moltypes {
            molecule_types.push(MoleculeType::parse(
                xdrfile,
                precision,
                tpr_version,
                symbol_table,
                ffparams,
            )?);
        }

        // get molecule blocks
        let n_molblocks = xdrfile.read_i32()?;

        let mut molecule_blocks = Vec::with_capacity(n_molblocks as usize);
        for _ in 0..n_molblocks {
            molecule_blocks.push(MolBlock::parse(xdrfile, precision)?)
        }

        // construct the topology from the molecule types and molecule blocks
        let topology = TprTopology::construct_topology(molecule_blocks, molecule_types)?;

        // read the number of atoms for sanity checking
        let n_atoms = xdrfile.read_i32()?;
        if n_atoms != expected_n_atoms {
            return Err(ParseTprError::InconsistentNumberOfAtoms(
                expected_n_atoms,
                n_atoms,
            ));
        }

        if n_atoms != topology.atoms.len() as i32 {
            return Err(ParseTprError::InconsistentNumberOfAtoms(
                expected_n_atoms,
                topology.atoms.len() as i32,
            ));
        }

        Ok(topology)
    }

    fn construct_topology(
        molecule_blocks: Vec<MolBlock>,
        molecule_types: Vec<MoleculeType>,
    ) -> Result<TprTopology, ParseTprError> {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        let mut atom_counter = 1;
        let mut residue_counter = 0;

        for molblock in molecule_blocks {
            let (new_atoms, new_bonds) = molblock.unpack2molecules(
                &molecule_types,
                &mut atom_counter,
                &mut residue_counter,
            )?;

            atoms.extend(new_atoms);
            bonds.extend(new_bonds);
        }

        Ok(TprTopology { atoms, bonds })
    }
}
