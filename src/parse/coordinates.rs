// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file contains functions for parsing positions, velocities and forces.

use crate::{errors::ParseTprError, Precision, TprCoordinates, TprHeader};

use super::xdr::XdrFile;

impl TprCoordinates {
    /// Get positions, velocities, and forces of particles from a tpr file.
    pub(super) fn parse(
        xdrfile: &mut XdrFile,
        tpr_header: &TprHeader,
    ) -> Result<TprCoordinates, ParseTprError> {
        let positions = if tpr_header.has_positions {
            TprCoordinates::read_block(xdrfile, tpr_header.precision, tpr_header.n_atoms)?
        } else {
            Vec::default()
        };

        let velocities = if tpr_header.has_velocities {
            TprCoordinates::read_block(xdrfile, tpr_header.precision, tpr_header.n_atoms)?
        } else {
            Vec::default()
        };

        let forces = if tpr_header.has_forces {
            TprCoordinates::read_block(xdrfile, tpr_header.precision, tpr_header.n_atoms)?
        } else {
            Vec::default()
        };

        Ok(TprCoordinates {
            positions,
            velocities,
            forces,
        })
    }

    /// Read a block of coordinates.
    fn read_block(
        xdrfile: &mut XdrFile,
        precision: Precision,
        n_items: i32,
    ) -> Result<Vec<[f64; 3]>, ParseTprError> {
        (0..n_items)
            .map(|_| xdrfile.read_vector3(precision))
            .collect::<Result<Vec<[f64; 3]>, std::io::Error>>()
            .map_err(ParseTprError::CouldNotRead)
    }
}
