// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file contains functions for parsing a tpr file.

use crate::{
    errors::ParseTprError,
    structures::{Precision, SimBox, TprFile, TprHeader, TprTopology},
    TprCoordinates,
};
use std::{fs::File, io::BufReader, path::Path};
use xdr::XdrFile;

use self::{ffparams::FFParams, symtab::SymTable};

pub mod coordinates;
pub mod ffparams;
pub mod header;
pub mod interactions;
pub mod molblocks;
pub mod moltypes;
pub mod simbox;
pub mod symtab;
pub mod topology;
pub mod xdr;

/// Parse a file in a Gromacs TPR format.
pub(crate) fn parse_tpr(filename: impl AsRef<Path>) -> Result<TprFile, ParseTprError> {
    let file = match File::open(filename.as_ref()) {
        Ok(x) => x,
        Err(_) => return Err(ParseTprError::CouldNotOpen(Box::from(filename.as_ref()))),
    };

    let reader = BufReader::new(file);
    let mut xdrfile = XdrFile::new(reader);

    // read header of the tpr file
    let header = TprHeader::parse(&mut xdrfile)?;

    // read simulation box (if present)
    let simbox = if header.has_box {
        Some(SimBox::parse(&mut xdrfile, header.precision)?)
    } else {
        None
    };

    // skip some data that used to be temperature coupling information
    let jump = match header.precision {
        Precision::Single => 4 * header.n_coupling_groups,
        Precision::Double => 8 * header.n_coupling_groups,
    };
    xdrfile.jump(jump as i64)?;

    // read symbol table
    let symtab = SymTable::parse(&mut xdrfile, header.tpr_version)?;

    // get system name
    let system_name = symtab.symstring(&mut xdrfile)?;

    // get force-field parameters
    let ffparams = FFParams::parse(&mut xdrfile, header.precision, header.tpr_version)?;

    let top = TprTopology::parse(
        &mut xdrfile,
        header.precision,
        header.tpr_version,
        &symtab,
        &ffparams,
        header.n_atoms,
    )?;

    // get positions, velocities, and forces
    let coordinates = TprCoordinates::parse(&mut xdrfile, &header)?;

    Ok(TprFile {
        header,
        system_name,
        simbox,
        topology: top,
        coordinates,
    })
}
