// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file contains functions for obtaining force-field parameters from TPR file.

use num;
use num_derive::FromPrimitive;
use std::collections::HashMap;
use strum::{EnumCount, EnumIter};

use crate::{errors::ParseTprError, structures::Precision};

use super::xdr::XdrFile;

/// Structure representing the force-field parameters of the TPR file.
/// `minitpr` ignores most of the parameters and only stores those that are important
/// for obtaining topology of the system.
#[derive(Debug, Clone)]
pub(super) struct FFParams {
    pub interaction_types: Vec<InteractionType>,
}

impl FFParams {
    /// Get `FFParams` from `XdrFile`.
    pub(super) fn parse(
        xdrfile: &mut XdrFile,
        precision: Precision,
        tpr_version: i32,
    ) -> Result<Self, ParseTprError> {
        // ignore the number of atom types
        xdrfile.jump(4)?;
        // get the number of function (interaction) types
        let n_interaction_types = xdrfile.read_i32()?;

        // get the interaction types
        let mut interaction_types = Vec::with_capacity(n_interaction_types as usize);
        for _ in 0..n_interaction_types {
            interaction_types.push(xdrfile.read_i32()?);
        }

        // reppow (always f64) and fudgeqq
        xdrfile.jump(8)?;
        xdrfile.skip_real(precision)?;

        let mut interaction_types_enum = Vec::with_capacity(n_interaction_types as usize);

        // renumber (update) the interaction types
        let updater = FTUpdater::default();

        // loop over all interaction types
        for interaction in interaction_types.iter_mut() {
            // loop over all file versions
            for (version, number) in &updater.update {
                // if the version of the tpr file is lower
                // and the function type is higher than the updated function type, renumber it
                if tpr_version < *version && *interaction >= *number {
                    *interaction += 1;
                }
            }

            let interaction_type_enum = num::FromPrimitive::from_i32(*interaction).expect(
                "FATAL MINITPR ERROR | FFParams::parse | Cannot convert interaction type to enum.",
            );

            interaction_types_enum.push(interaction_type_enum);

            // get parameters of the function type
            Self::get_params(xdrfile, interaction_type_enum, precision, tpr_version)?;
        }

        Ok(FFParams {
            interaction_types: interaction_types_enum,
        })
    }

    /// Read parameters for the target interaction type from the xdr file.
    /// This function does not return anything, if successful.
    /// The parameters are read and then promptly ignored as we do not need them.
    fn get_params(
        xdrfile: &mut XdrFile,
        interaction_type: InteractionType,
        precision: Precision,
        tpr_version: i32,
    ) -> Result<(), ParseTprError> {
        match interaction_type {
            InteractionType::F_ANGLES
            | InteractionType::F_G96ANGLES
            | InteractionType::F_BONDS
            | InteractionType::F_G96BONDS
            | InteractionType::F_HARMONIC
            | InteractionType::F_IDIHS => {
                xdrfile.skip_multiple_reals(precision, 4)?;
            }
            InteractionType::F_RESTRANGLES => {
                xdrfile.skip_multiple_reals(precision, 2)?;
                if tpr_version >= 134 {
                    xdrfile.skip_multiple_reals(precision, 2)?;
                }
            }
            InteractionType::F_LINEAR_ANGLES => {
                xdrfile.skip_multiple_reals(precision, 4)?;
            }
            InteractionType::F_FENEBONDS => {
                xdrfile.skip_multiple_reals(precision, 2)?;
            }
            InteractionType::F_RESTRBONDS => {
                xdrfile.skip_multiple_reals(precision, 8)?;
            }
            InteractionType::F_TABBONDS
            | InteractionType::F_TABBONDSNC
            | InteractionType::F_TABANGLES
            | InteractionType::F_TABDIHS => {
                xdrfile.skip_multiple_reals(precision, 2)?;
                xdrfile.jump(4)?;
            }
            InteractionType::F_CROSS_BOND_BONDS => {
                xdrfile.skip_multiple_reals(precision, 3)?;
            }
            InteractionType::F_CROSS_BOND_ANGLES => {
                xdrfile.skip_multiple_reals(precision, 4)?;
            }
            InteractionType::F_UREY_BRADLEY => {
                xdrfile.skip_multiple_reals(precision, 8)?;
            }
            InteractionType::F_QUARTIC_ANGLES => {
                xdrfile.skip_multiple_reals(precision, 6)?;
            }
            InteractionType::F_BHAM => {
                xdrfile.skip_multiple_reals(precision, 3)?;
            }
            InteractionType::F_MORSE => {
                xdrfile.skip_multiple_reals(precision, 6)?;
            }
            InteractionType::F_CUBICBONDS => {
                xdrfile.skip_multiple_reals(precision, 3)?;
            }
            InteractionType::F_CONNBONDS => {}
            InteractionType::F_POLARIZATION => {
                xdrfile.skip_real(precision)?;
            }
            InteractionType::F_ANHARM_POL => {
                xdrfile.skip_multiple_reals(precision, 3)?;
            }
            InteractionType::F_WATER_POL => {
                xdrfile.skip_multiple_reals(precision, 6)?;
            }
            InteractionType::F_THOLE_POL => {
                xdrfile.skip_multiple_reals(precision, 3)?;
                if tpr_version < 127 {
                    xdrfile.skip_real(precision)?;
                }
            }
            InteractionType::F_LJ => {
                xdrfile.skip_multiple_reals(precision, 2)?;
            }
            InteractionType::F_LJ14 => {
                xdrfile.skip_multiple_reals(precision, 4)?;
            }
            InteractionType::F_LJC14_Q => {
                xdrfile.skip_multiple_reals(precision, 5)?;
            }
            InteractionType::F_LJC_PAIRS_NB => {
                xdrfile.skip_multiple_reals(precision, 4)?;
            }
            InteractionType::F_PDIHS
            | InteractionType::F_PIDIHS
            | InteractionType::F_ANGRES
            | InteractionType::F_ANGRESZ => {
                xdrfile.skip_multiple_reals(precision, 4)?;
                xdrfile.jump(4)?;
            }
            InteractionType::F_RESTRDIHS => {
                xdrfile.skip_multiple_reals(precision, 2)?;
                if tpr_version >= 134 {
                    xdrfile.skip_multiple_reals(precision, 2)?;
                }
            }
            InteractionType::F_DISRES => {
                xdrfile.jump(8)?;
                xdrfile.skip_multiple_reals(precision, 4)?;
            }
            InteractionType::F_ORIRES => {
                xdrfile.jump(12)?;
                xdrfile.skip_multiple_reals(precision, 3)?;
            }
            InteractionType::F_DIHRES => {
                xdrfile.skip_multiple_reals(precision, 6)?;
            }
            InteractionType::F_POSRES => {
                xdrfile.skip_multiple_reals(precision, 4 * crate::DIM as i64)?;
            }
            InteractionType::F_FBPOSRES => {
                xdrfile.jump(4)?;
                xdrfile.skip_multiple_reals(precision, 2 + crate::DIM as i64)?;
            }
            InteractionType::F_CBTDIHS => {
                xdrfile.skip_multiple_reals(precision, crate::NR_CBTDIHS as i64)?;
                if tpr_version >= 134 {
                    xdrfile.skip_multiple_reals(precision, crate::NR_CBTDIHS as i64)?;
                }
            }
            InteractionType::F_RBDIHS | InteractionType::F_FOURDIHS => {
                xdrfile.skip_multiple_reals(precision, 2 * crate::NR_RBDIHS as i64)?;
            }
            InteractionType::F_CONSTR | InteractionType::F_CONSTRNC => {
                xdrfile.skip_multiple_reals(precision, 2)?;
            }
            InteractionType::F_SETTLE => {
                xdrfile.skip_multiple_reals(precision, 2)?;
            }
            InteractionType::F_VSITE1 => {}
            InteractionType::F_VSITE2 | InteractionType::F_VSITE2FD => {
                xdrfile.skip_real(precision)?;
            }
            InteractionType::F_VSITE3
            | InteractionType::F_VSITE3FD
            | InteractionType::F_VSITE3FAD => {
                xdrfile.skip_multiple_reals(precision, 2)?;
            }
            InteractionType::F_VSITE3OUT
            | InteractionType::F_VSITE4FD
            | InteractionType::F_VSITE4FDN => {
                xdrfile.skip_multiple_reals(precision, 3)?;
            }
            InteractionType::F_VSITEN => {
                xdrfile.jump(4)?;
                xdrfile.skip_real(precision)?;
            }
            InteractionType::F_GB12_NOLONGERUSED
            | InteractionType::F_GB13_NOLONGERUSED
            | InteractionType::F_GB14_NOLONGERUSED => {
                if tpr_version < 113 {
                    xdrfile.skip_multiple_reals(precision, 5)?;
                }
            }
            InteractionType::F_CMAP => {
                xdrfile.jump(8)?;
            }
            // Ignore the other function types...
            _ => (),
        }

        Ok(())
    }
}

/// Takes care of updating the interaction type numbers.
pub(super) struct FTUpdater {
    /// format is `file version`: `interaction type number`
    pub update: HashMap<i32, i32>,
}

impl FTUpdater {
    /// Create a `TprFTUpdater` structure.
    pub(super) fn default() -> Self {
        FTUpdater {
            update: HashMap::from([(121, 65), (118, 67), (117, 76), (137, 78)]),
        }
    }
}

/// Enum describing all supported interaction types.
#[derive(Debug, Clone, Copy, FromPrimitive, EnumIter, EnumCount)]
#[allow(non_camel_case_types, dead_code)]
pub(crate) enum InteractionType {
    F_BONDS = 0,
    F_G96BONDS,
    F_MORSE,
    F_CUBICBONDS,
    F_CONNBONDS,
    F_HARMONIC,
    F_FENEBONDS,
    F_TABBONDS,
    F_TABBONDSNC,
    F_RESTRBONDS,
    F_ANGLES,
    F_G96ANGLES,
    F_RESTRANGLES,
    F_LINEAR_ANGLES,
    F_CROSS_BOND_BONDS,
    F_CROSS_BOND_ANGLES,
    F_UREY_BRADLEY,
    F_QUARTIC_ANGLES,
    F_TABANGLES,
    F_PDIHS,
    F_RBDIHS,
    F_RESTRDIHS,
    F_CBTDIHS,
    F_FOURDIHS,
    F_IDIHS,
    F_PIDIHS,
    F_TABDIHS,
    F_CMAP,
    F_GB12_NOLONGERUSED,
    F_GB13_NOLONGERUSED,
    F_GB14_NOLONGERUSED,
    F_GBPOL_NOLONGERUSED,
    F_NPSOLVATION_NOLONGERUSED,
    F_LJ14,
    F_COUL14,
    F_LJC14_Q,
    F_LJC_PAIRS_NB,
    F_LJ,
    F_BHAM,
    F_LJ_LR_NOLONGERUSED,
    F_BHAM_LR_NOLONGERUSED,
    F_DISPCORR,
    F_COUL_SR,
    F_COUL_LR_NOLONGERUSED,
    F_RF_EXCL,
    F_COUL_RECIP,
    F_LJ_RECIP,
    F_DPD,
    F_POLARIZATION,
    F_WATER_POL,
    F_THOLE_POL,
    F_ANHARM_POL,
    F_POSRES,
    F_FBPOSRES,
    F_DISRES,
    F_DISRESVIOL,
    F_ORIRES,
    F_ORIRESDEV,
    F_ANGRES,
    F_ANGRESZ,
    F_DIHRES,
    F_DIHRESVIOL,
    F_CONSTR,
    F_CONSTRNC,
    F_SETTLE,
    F_VSITE1,
    F_VSITE2,
    F_VSITE2FD,
    F_VSITE3,
    F_VSITE3FD,
    F_VSITE3FAD,
    F_VSITE3OUT,
    F_VSITE4FD,
    F_VSITE4FDN,
    F_VSITEN,
    F_COM_PULL,
    F_DENSITYFITTING,
    F_EQM,
    F_ENNPOT,
    F_EPOT,
    F_EKIN,
    F_ETOT,
    F_ECONSERVED,
    F_TEMP,
    F_VTEMP_NOLONGERUSED,
    F_PDISPCORR,
    F_PRES,
    F_DVDL_CONSTR,
    F_DVDL,
    F_DKDL,
    F_DVDL_COUL,
    F_DVDL_VDW,
    F_DVDL_BONDED,
    F_DVDL_RESTRAINT,
    F_DVDL_TEMPERATURE,
}

impl InteractionType {
    /// Get the number of interacting atoms for this InteractionType.
    pub(super) fn n_interacting_atoms(&self) -> i32 {
        match self {
            InteractionType::F_POSRES | InteractionType::F_FBPOSRES => 1,

            InteractionType::F_BONDS
            | InteractionType::F_G96BONDS
            | InteractionType::F_MORSE
            | InteractionType::F_CUBICBONDS
            | InteractionType::F_CONNBONDS
            | InteractionType::F_HARMONIC
            | InteractionType::F_FENEBONDS
            | InteractionType::F_TABBONDS
            | InteractionType::F_TABBONDSNC
            | InteractionType::F_RESTRBONDS
            | InteractionType::F_GB12_NOLONGERUSED
            | InteractionType::F_GB13_NOLONGERUSED
            | InteractionType::F_GB14_NOLONGERUSED
            | InteractionType::F_LJ14
            | InteractionType::F_LJC14_Q
            | InteractionType::F_LJC_PAIRS_NB
            | InteractionType::F_LJ
            | InteractionType::F_BHAM
            | InteractionType::F_POLARIZATION
            | InteractionType::F_ANHARM_POL
            | InteractionType::F_DISRES
            | InteractionType::F_ORIRES
            | InteractionType::F_ANGRESZ
            | InteractionType::F_CONSTR
            | InteractionType::F_CONSTRNC
            | InteractionType::F_VSITE1
            | InteractionType::F_VSITEN => 2,

            InteractionType::F_ANGLES
            | InteractionType::F_G96ANGLES
            | InteractionType::F_RESTRANGLES
            | InteractionType::F_LINEAR_ANGLES
            | InteractionType::F_CROSS_BOND_BONDS
            | InteractionType::F_CROSS_BOND_ANGLES
            | InteractionType::F_UREY_BRADLEY
            | InteractionType::F_QUARTIC_ANGLES
            | InteractionType::F_TABANGLES
            | InteractionType::F_SETTLE
            | InteractionType::F_VSITE2
            | InteractionType::F_VSITE2FD => 3,

            InteractionType::F_PDIHS
            | InteractionType::F_RBDIHS
            | InteractionType::F_RESTRDIHS
            | InteractionType::F_CBTDIHS
            | InteractionType::F_FOURDIHS
            | InteractionType::F_IDIHS
            | InteractionType::F_PIDIHS
            | InteractionType::F_TABDIHS
            | InteractionType::F_THOLE_POL
            | InteractionType::F_ANGRES
            | InteractionType::F_DIHRES
            | InteractionType::F_VSITE3
            | InteractionType::F_VSITE3FD
            | InteractionType::F_VSITE3FAD
            | InteractionType::F_VSITE3OUT => 4,

            InteractionType::F_CMAP
            | InteractionType::F_WATER_POL
            | InteractionType::F_VSITE4FD
            | InteractionType::F_VSITE4FDN => 5,

            _ => 0,
        }
    }
}
