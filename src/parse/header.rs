// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file contains functions for parsing the TPR file header.

use crate::{
    errors::ParseTprError,
    structures::{Precision, TprHeader},
};

use super::xdr::XdrFile;

impl TprHeader {
    /// Get TprHeader from a tpr file.
    pub(super) fn parse(xdrfile: &mut XdrFile) -> Result<TprHeader, ParseTprError> {
        // get gromacs version used to write the tpr file
        let gromacs_version = xdrfile.read_string_4byte()?;

        // check that this is indeed a tpr file
        if !gromacs_version.contains("VERSION") {
            return Err(ParseTprError::NotTpr);
        }

        // get precision of the tpr file
        let precision = match xdrfile.read_i32()? {
            4 => Precision::Single,
            8 => Precision::Double,
            x => return Err(ParseTprError::UnsupportedPrecision(x)),
        };

        // get version of the file
        let tpr_version = xdrfile.read_i32()?;

        // check that the version of the tpr file is supported
        if tpr_version < 103 {
            return Err(ParseTprError::UnsupportedVersion(tpr_version));
        }

        let tpr_generation = xdrfile.read_i32()?;
        let file_tag = xdrfile.read_string_4byte()?;
        let n_atoms = xdrfile.read_i32()?;
        let n_coupling_groups = xdrfile.read_i32()?;
        let fep_state = xdrfile.read_i32()?;
        let lambda = xdrfile.read_real(precision)?;

        let has_input_record = xdrfile.read_bool_header()?;
        let has_topology = xdrfile.read_bool_header()?;
        let has_coordinates = xdrfile.read_bool_header()?;
        let has_velocities = xdrfile.read_bool_header()?;
        let has_forces = xdrfile.read_bool_header()?;
        let has_box = xdrfile.read_bool_header()?;

        // get information about the size of the tpr file (version >= 119 && generation >= 27)
        let body_size = if tpr_version >= 119 && tpr_generation >= 27 {
            Some(xdrfile.read_i64()?)
        } else {
            None
        };

        Ok(TprHeader {
            gromacs_version,
            precision,
            tpr_version,
            tpr_generation,
            file_tag,
            n_atoms,
            n_coupling_groups,
            fep_state,
            lambda,
            has_input_record,
            has_topology,
            has_coordinates,
            has_velocities,
            has_forces,
            has_box,
            body_size,
        })
    }
}
