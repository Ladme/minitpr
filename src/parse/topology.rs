// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file contains functions for obtaining system topology from a TPR file.

use super::{
    coordinates::Coordinates, ffparams::FFParams, interactions::Interaction, molblocks::MolBlock,
    moltypes::MoleculeType, xdr::XdrFile,
};
use crate::{
    errors::ParseTprError,
    structures::{Precision, TprTopology},
    NR_GROUP_TYPES,
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
    ) -> Result<Self, ParseTprError> {
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

        // read the number of atoms for sanity checking
        let n_atoms = xdrfile.read_i32()?;

        // read intermolecular interactions
        let intermolecular = if xdrfile.read_bool_body(tpr_version)? {
            Some(super::interactions::read_interactions(
                xdrfile,
                tpr_version,
                ffparams,
            )?)
        } else {
            None
        };

        // construct the topology from the molecule types, molecule blocks and intermolecular interactions
        let topology =
            TprTopology::construct_topology(molecule_blocks, molecule_types, intermolecular)?;

        // check that the number of atoms is consistent
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

        // skip atom types
        if tpr_version < 128 {
            let n_types = xdrfile.read_i32()?;
            if tpr_version < 113 {
                xdrfile.skip_multiple_reals(precision, 5 * n_types as i64)?;
            }

            xdrfile.jump(4 * n_types as i64)?;
        }

        // skip dihedral correction maps
        let n_grids = xdrfile.read_i32()?;
        let grid_spacing = xdrfile.read_i32()?;
        xdrfile.skip_multiple_reals(
            precision,
            (4 * n_grids * grid_spacing * grid_spacing) as i64,
        )?;

        // skip atom groups
        for _ in 0..NR_GROUP_TYPES {
            let group_size = xdrfile.read_i32()?;
            xdrfile.jump(4 * group_size as i64)?;
        }

        let n_group_names = xdrfile.read_i32()?;
        xdrfile.jump(4 * n_group_names as i64)?;

        for _ in 0..NR_GROUP_TYPES {
            let n_group_numbers = xdrfile.read_i32()?;
            xdrfile.skip_multiple_uchars_body(tpr_version, n_group_numbers as i64)?;
        }

        // skip exclusions
        if tpr_version >= 120 {
            let intermolecular_exclusion_group_size = xdrfile.read_i64()?;
            if intermolecular_exclusion_group_size < 0 {
                return Err(ParseTprError::InvalidIntermolecularExclusionGroupSize(
                    intermolecular_exclusion_group_size,
                ));
            }

            xdrfile.jump(4 * intermolecular_exclusion_group_size)?;
        }

        Ok(topology)
    }

    /// Construct the final topology from molecule blocks, molecule types and intermolecular interactions.
    fn construct_topology(
        molecule_blocks: Vec<MolBlock>,
        molecule_types: Vec<MoleculeType>,
        intermolecular: Option<Vec<Interaction>>,
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

        // convert intermolecular interactions to bonds
        if let Some(inter) = intermolecular {
            for interaction in inter.iter() {
                if let Some(bond) = interaction.unpack2bond(&atoms)? {
                    bonds.push(bond);
                }
            }
        }

        Ok(TprTopology { atoms, bonds })
    }

    /// Get positions, velocities, and forces for particles in the topology from the `Coordinates` structure.
    pub(super) fn fill_with_coordinates(&mut self, coordinates: Coordinates) {
        for (pos, atom) in coordinates.positions.into_iter().zip(self.atoms.iter_mut()) {
            atom.position = Some(pos);
        }

        for (vel, atom) in coordinates
            .velocities
            .into_iter()
            .zip(self.atoms.iter_mut())
        {
            atom.velocity = Some(vel);
        }

        for (force, atom) in coordinates.forces.into_iter().zip(self.atoms.iter_mut()) {
            atom.force = Some(force);
        }
    }
}
