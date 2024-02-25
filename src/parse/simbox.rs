// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file contains functions for parsing simulation box information from a TPR file.

use crate::structures::Precision;
use crate::DIM;
use crate::{errors::ParseTprError, structures::SimBox};

use super::xdr::XdrFile;

impl SimBox {
    /// Get simulation box dimensions of specified precision from the tpr file.
    pub(super) fn parse(
        xdrfile: &mut XdrFile,
        precision: Precision,
    ) -> Result<SimBox, ParseTprError> {
        fn fill_matrix(
            xdrfile: &mut XdrFile,
            precision: Precision,
        ) -> Result<[[f64; DIM]; DIM], ParseTprError> {
            let mut matrix = [[0.0f64; DIM]; DIM];
            for fields in matrix.iter_mut().take(DIM) {
                for field in fields.iter_mut().take(DIM) {
                    *field = xdrfile.read_real(precision)?;
                }
            }
            Ok(matrix)
        }

        let simbox = fill_matrix(xdrfile, precision)?;
        let simbox_rel = fill_matrix(xdrfile, precision)?;
        let simbox_v = fill_matrix(xdrfile, precision)?;

        Ok(SimBox {
            simbox,
            simbox_rel,
            simbox_v,
        })
    }
}
