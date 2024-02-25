// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file contains functions for working with Symbol Table.

use super::xdr::XdrFile;
use crate::errors::ParseTprError;

/// Structure representing the Symbol Table.
#[derive(Debug, Clone)]
pub(super) struct SymTable {
    pub symbols: Vec<String>,
}

impl SymTable {
    /// Get `SymTable` from `XdrFile`.
    pub(super) fn parse(xdrfile: &mut XdrFile, tpr_version: i32) -> Result<Self, ParseTprError> {
        let symtab_len = xdrfile.read_i32()?;

        let mut symtab = SymTable {
            symbols: Vec::with_capacity(symtab_len as usize),
        };

        for _ in 0..symtab_len {
            symtab.symbols.push(xdrfile.read_string_body(tpr_version)?);
        }

        Ok(symtab)
    }

    /// Read `i32` from `XdrFile` and convert it to string using the `SymTable`.
    pub(super) fn symstring(&self, xdrfile: &mut XdrFile) -> Result<String, ParseTprError> {
        let index = xdrfile.read_i32()?;

        Ok(match self.symbols.get(index as usize) {
            Some(x) => x,
            None => return Err(ParseTprError::IndexNotInSymTable(index)),
        }
        .to_owned())
    }
}
