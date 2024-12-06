// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file contains functions for parsing intramolecular and intermolecular interactions.

use strum::IntoEnumIterator;

use crate::{errors::ParseTprError, Atom, Bond};

use super::{
    ffparams::{FFParams, FTUpdater, InteractionType},
    xdr::XdrFile,
};

/// Structure representing a intramolecular or an intermolecular interaction.
#[derive(Debug, Clone)]
pub(super) struct Interaction {
    pub interaction_type: InteractionType,
    pub interacting_atom_indices: Vec<i32>,
}

/// Read intramolecular or intermolecular interactions.
pub(super) fn read_interactions(
    xdrfile: &mut XdrFile,
    tpr_version: i32,
    ffparams: &FFParams,
) -> Result<Vec<Interaction>, ParseTprError> {
    let updater = FTUpdater::default();
    let mut interactions = Vec::new();

    for functype in InteractionType::iter() {
        let mut read = true;

        for (version, number) in &updater.update {
            if tpr_version < *version && functype as i32 == *number {
                read = false;
                break;
            }
        }

        if !read {
            continue;
        }

        let number_of_instances = xdrfile.read_i32()?;
        // get the number of atoms interacting via this interaction type
        let n_interacting_atoms = functype.n_interacting_atoms();

        // sanity check: the number of instances must be divisible by the number of atoms + 1
        if number_of_instances % (n_interacting_atoms + 1) != 0 {
            return Err(ParseTprError::InteractionDiscrepancy(functype as i32));
        }

        for _ in (0..number_of_instances).step_by(n_interacting_atoms as usize + 1) {
            interactions.push(Interaction::parse(xdrfile, n_interacting_atoms, ffparams)?);
        }
    }

    Ok(interactions)
}

impl Interaction {
    /// Get `Interaction` from an `XdrFile`.
    fn parse(
        xdrfile: &mut XdrFile,
        n_interacting_atoms: i32,
        ffparams: &FFParams,
    ) -> Result<Self, ParseTprError> {
        let interaction_type_index = xdrfile.read_i32()?;
        let interaction_type = match ffparams
            .interaction_types
            .get(interaction_type_index as usize)
        {
            Some(x) => *x,
            None => {
                return Err(ParseTprError::InvalidInteractionType(
                    interaction_type_index,
                ))
            }
        };

        let mut interacting_atom_indices = Vec::with_capacity(n_interacting_atoms as usize);

        for _ in 0..n_interacting_atoms {
            interacting_atom_indices.push(xdrfile.read_i32()?);
        }

        Ok(Interaction {
            interaction_type,
            interacting_atom_indices,
        })
    }

    /// Return `true` if the `Interaction` is considered to be a bond.
    /// Otherwise, return `false`.
    pub(super) fn is_bond(&self) -> bool {
        matches!(
            self.interaction_type,
            InteractionType::F_BONDS
                | InteractionType::F_G96BONDS
                | InteractionType::F_MORSE
                | InteractionType::F_CUBICBONDS
                | InteractionType::F_CONNBONDS
                | InteractionType::F_HARMONIC
                | InteractionType::F_FENEBONDS
                | InteractionType::F_RESTRBONDS
                | InteractionType::F_CONSTR
                | InteractionType::F_CONSTRNC
                | InteractionType::F_TABBONDS
                | InteractionType::F_TABBONDSNC
        )
    }

    /// Unpack SETTLE interaction into bonds.
    /// Returns an empty vector, if the interaction is not a settle.
    /// Returns `ParseTprError` if the bonds could not be constructed due to some inconsistency in the input data.
    pub(super) fn settle2bonds(&self, atoms: &[Atom]) -> Result<Vec<Bond>, ParseTprError> {
        if !matches!(self.interaction_type, InteractionType::F_SETTLE) {
            return Ok(vec![]);
        }

        // three atoms must be involved
        if self.interacting_atom_indices.len() != 3 {
            return Err(ParseTprError::InvalidNumberOfSettleAtoms(
                self.interacting_atom_indices.len(),
            ));
        }

        // get global atom indices
        let get_atom_index = |index: usize| -> Result<usize, ParseTprError> {
            atoms
                .get(self.interacting_atom_indices[index] as usize)
                .map(|x| (x.atom_number - 1) as usize)
                .ok_or(ParseTprError::CouldNotConstructTopology)
        };

        Ok(vec![
            Bond {
                atom1: get_atom_index(0)?,
                atom2: get_atom_index(1)?,
            },
            Bond {
                atom1: get_atom_index(0)?,
                atom2: get_atom_index(2)?,
            },
        ])
    }

    /// Unpack `Interaction` into an Bond between specific atoms.
    /// Returns `None`, if the interaction is not a bond.
    /// Returns `ParseTprError` if the Bond could not be constructed due to some inconsistency in the input data.
    pub(super) fn unpack2bond(&self, atoms: &[Atom]) -> Result<Option<Bond>, ParseTprError> {
        // check whether this interaction is a bond
        if !self.is_bond() {
            return Ok(None);
        }

        // bond must involve exactly two atoms
        if self.interacting_atom_indices.len() != 2 {
            return Err(ParseTprError::InvalidNumberOfBondedAtoms(
                self.interacting_atom_indices.len(),
            ));
        }

        // get global atom indices
        let get_atom_index = |index: usize| -> Result<usize, ParseTprError> {
            atoms
                .get(self.interacting_atom_indices[index] as usize)
                .map(|x| (x.atom_number - 1) as usize)
                .ok_or(ParseTprError::CouldNotConstructTopology)
        };

        Ok(Some(Bond {
            atom1: get_atom_index(0)?,
            atom2: get_atom_index(1)?,
        }))
    }
}
