// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file contains functions for obtaining molecule types from TPR file.

use mendeleev::Element;

use crate::{
    errors::ParseTprError,
    parse::xdr::XdrFile,
    structures::{Atom, Bond, Precision},
};

use super::{
    ffparams::FFParams,
    interactions::{self, Interaction},
    symtab::SymTable,
};

/// Structure representing Molecule Type.
#[derive(Debug, Clone)]
pub(super) struct MoleculeType {
    pub atoms: Vec<MoleculeTypeAtom>,
    pub residues: Vec<MoleculeTypeResidue>,
    pub interactions: Vec<Interaction>,
}

/// Structure representing an atom of a Molecule Type.
#[derive(Debug, Clone)]
pub(super) struct MoleculeTypeAtom {
    pub name: String,
    pub mass: f64,
    pub charge: f64,
    pub residue_index: i32,
    pub element: Option<Element>,
}

/// Structure representing a residue of a Molecule Type.
#[derive(Debug, Clone)]
pub(super) struct MoleculeTypeResidue {
    pub name: String,
    pub number: i32,
}

impl MoleculeType {
    /// Get `MoleculeType` from `XdrFile`.
    pub(super) fn parse(
        xdrfile: &mut XdrFile,
        precision: Precision,
        tpr_version: i32,
        symbol_table: &SymTable,
        ffparams: &FFParams,
    ) -> Result<Self, ParseTprError> {
        // skip the name of the molecule type
        symbol_table.symstring(xdrfile)?;

        // get the number of atoms and residues in the molecule type
        let n_atoms = xdrfile.read_i32()?;
        let n_residues = xdrfile.read_i32()?;

        // read atoms
        let mut atoms = Vec::with_capacity(n_atoms as usize);
        for _ in 0..(n_atoms as usize) {
            atoms.push(MoleculeTypeAtom::parse(xdrfile, precision, tpr_version)?);
        }

        // read atom names
        for atom in atoms.iter_mut() {
            atom.name = symbol_table.symstring(xdrfile)?;
        }

        // skip names and B names of the atom types
        for _ in atoms.iter() {
            symbol_table.symstring(xdrfile)?;
            symbol_table.symstring(xdrfile)?;
        }

        // read residues
        let mut residues = Vec::with_capacity(n_residues as usize);
        for _ in 0..(n_residues as usize) {
            residues.push(MoleculeTypeResidue::parse(
                xdrfile,
                tpr_version,
                symbol_table,
            )?);
        }

        // read interactions
        let interactions = interactions::read_interactions(xdrfile, tpr_version, ffparams)?;

        // skip block indices
        let n_blocks = xdrfile.read_i32()?;
        xdrfile.jump(4 * (n_blocks as i64 + 1))?;

        // skip exclusions
        let n_exclusions = xdrfile.read_i32()?;
        let n_excluded = xdrfile.read_i32()?;
        xdrfile.jump(4 * n_exclusions as i64 + 4)?;
        xdrfile.jump(4 * n_excluded as i64)?;

        Ok(MoleculeType {
            atoms,
            residues,
            interactions,
        })
    }

    /// Unpack `MoleculeType` to molecule, i.e., a vector of atoms and a vector of bonds.
    pub(super) fn unpack2molecule(
        &self,
        atom_counter: &mut i32,
        residue_counter: &mut i32,
    ) -> Result<(Vec<Atom>, Vec<Bond>), ParseTprError> {
        let mut atoms = Vec::with_capacity(self.atoms.len());

        let mut previous_residue_number = None;
        for moltype_atom in &self.atoms {
            atoms.push(moltype_atom.convert2atom(
                &self.residues,
                atom_counter,
                residue_counter,
                &mut previous_residue_number,
            )?)
        }

        let mut bonds = Vec::new();
        for interaction in self.interactions.iter() {
            match interaction.unpack2bond(&atoms) {
                Ok(Some(x)) => bonds.push(x),
                Ok(None) => (),
                Err(e) => return Err(e),
            }
        }

        Ok((atoms, bonds))
    }
}

impl MoleculeTypeAtom {
    /// Get `MoleculeTypeAtom` from an `XdrFile`.
    fn parse(
        xdrfile: &mut XdrFile,
        precision: Precision,
        tpr_version: i32,
    ) -> Result<Self, ParseTprError> {
        let mass = xdrfile.read_real(precision)?;
        let charge = xdrfile.read_real(precision)?;

        // ignore mass_b and charge_b
        xdrfile.skip_multiple_reals(precision, 2)?;

        // skip both atom type indices
        xdrfile.read_ushort_body(tpr_version)?;
        xdrfile.read_ushort_body(tpr_version)?;

        // skip p-type
        xdrfile.jump(4)?;
        let residue_index = xdrfile.read_i32()?;

        let atomic_number = xdrfile.read_i32()?;
        let element = from_atom_number(atomic_number);

        Ok(MoleculeTypeAtom {
            name: String::from("Unknown"),
            mass,
            charge,
            residue_index,
            element,
        })
    }

    /// Convert `MoleculeTypeAtom` to `Atom` structure.
    pub fn convert2atom(
        &self,
        residues: &[MoleculeTypeResidue],
        atom_counter: &mut i32,
        residue_counter: &mut i32,
        previous_residue_number: &mut Option<i32>,
    ) -> Result<Atom, ParseTprError> {
        let residue = match residues.get(self.residue_index as usize) {
            Some(x) => x,
            None => return Err(ParseTprError::CouldNotConstructTopology),
        };

        // increase the residue counter, if new residue is encountered
        if let Some(unwrapped) = previous_residue_number {
            if *unwrapped != residue.number {
                *residue_counter += 1;
                *previous_residue_number = Some(residue.number);
            }
        } else {
            *residue_counter += 1;
            *previous_residue_number = Some(residue.number);
        }

        *atom_counter += 1;

        Ok(Atom {
            atom_name: self.name.clone(),
            atom_number: *atom_counter - 1,
            residue_name: residue.name.clone(),
            residue_number: *residue_counter,
            mass: self.mass,
            charge: self.charge,
            element: self.element,
            position: None,
            velocity: None,
            force: None,
        })
    }
}

impl MoleculeTypeResidue {
    /// Get `MoleculeTypeResidue` from an `XdrFile`.
    fn parse(
        xdrfile: &mut XdrFile,
        tpr_version: i32,
        symbol_table: &SymTable,
    ) -> Result<Self, ParseTprError> {
        let name = symbol_table.symstring(xdrfile)?;
        let number = xdrfile.read_i32()?;

        // skip insertion code
        xdrfile.read_uchar_body(tpr_version)?;

        Ok(MoleculeTypeResidue { name, number })
    }
}

/// Attempt to identify element from atomic number.
/// Returns `None` if the atomic number corresponds to no element.
fn from_atom_number(atomic_number: i32) -> Option<Element> {
    Element::list().get(atomic_number as usize - 1).copied()
}
