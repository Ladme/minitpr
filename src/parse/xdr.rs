// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

//! This file contains low-level functions for reading XDR files.

use std::{
    fs::File,
    io::{BufReader, Error, Read},
};

use byteorder::{BigEndian, ReadBytesExt};

use crate::structures::Precision;

/// Structure representing the TPR file being read.
#[derive(Debug)]
pub(super) struct XdrFile {
    reader: BufReader<File>,
}

impl XdrFile {
    /// Create a new `XdrFile` structure.
    #[inline(always)]
    pub(super) fn new(reader: BufReader<File>) -> Self {
        XdrFile { reader }
    }

    /// Jump forward by N bytes.
    #[inline(always)]
    pub(super) fn jump(&mut self, n_bytes: i64) -> Result<(), Error> {
        self.reader.seek_relative(n_bytes)
    }

    /// Read `u8` value from `XdrFile`.
    #[inline(always)]
    pub(super) fn read_u8(&mut self) -> Result<u8, Error> {
        self.reader.read_u8()
    }

    /// Read `u16` value from `XdrFile`.
    #[inline(always)]
    pub(super) fn read_u16(&mut self) -> Result<u16, Error> {
        self.reader.read_u16::<BigEndian>()
    }

    /// Read `i32` value from `XdrFile`.
    #[inline(always)]
    pub(super) fn read_i32(&mut self) -> Result<i32, Error> {
        self.reader.read_i32::<BigEndian>()
    }

    /// Read `u32` value from `XdrFile`.
    #[inline(always)]
    pub(super) fn read_u32(&mut self) -> Result<u32, Error> {
        self.reader.read_u32::<BigEndian>()
    }

    /// Read `u64` value from `XdrFile`.
    #[inline(always)]
    pub(super) fn read_u64(&mut self) -> Result<u64, Error> {
        self.reader.read_u64::<BigEndian>()
    }

    /// Read `i64` value from `XdrFile`.
    #[inline(always)]
    pub(super) fn read_i64(&mut self) -> Result<i64, Error> {
        self.reader.read_i64::<BigEndian>()
    }

    /// Read `f32` value from `XdrFile`.
    #[inline(always)]
    pub(super) fn read_f32(&mut self) -> Result<f32, Error> {
        self.reader.read_f32::<BigEndian>()
    }

    /// Read `f64` value from `XdrFile`.
    #[inline(always)]
    pub(super) fn read_f64(&mut self) -> Result<f64, Error> {
        self.reader.read_f64::<BigEndian>()
    }

    /// Read `f32` value or `f64` value from `XdrFile` depending on the provided precision.
    #[inline(always)]
    pub(super) fn read_real(&mut self, precision: Precision) -> Result<f64, Error> {
        match precision {
            Precision::Single => Ok(self.read_f32()? as f64),
            Precision::Double => self.read_f64(),
        }
    }

    /// Jump N bytes depending on the provided precision.
    #[inline(always)]
    pub(super) fn skip_real(&mut self, precision: Precision) -> Result<(), Error> {
        match precision {
            Precision::Single => self.jump(4),
            Precision::Double => self.jump(8),
        }
    }

    /// Jump N bytes depending on the provided precision and the number of real numbers to skip.
    #[inline(always)]
    pub(super) fn skip_multiple_reals(
        &mut self,
        precision: Precision,
        n_reals: i64,
    ) -> Result<(), Error> {
        match precision {
            Precision::Single => self.jump(4 * n_reals),
            Precision::Double => self.jump(8 * n_reals),
        }
    }

    /// Read `bool` value from `XdrFile` as an `u32` value. This function is used ONLY in the TPR header.
    #[inline(always)]
    pub(super) fn read_bool_header(&mut self) -> Result<bool, Error> {
        Ok(self.read_u32()? != 0)
    }

    /// Read a string with one useless 4byte header and one useful 4byte header from `XdrFile`.
    /// This is used for a) the tpr file header and b) for the body of tpr files version < 119.
    pub(super) fn read_string_4byte(&mut self) -> Result<String, Error> {
        // first 4 bytes of the string header are not used
        self.reader.seek_relative(4)?;

        // get length of the string
        let mut len = self.read_u32()?;

        // make sure that the length is a multiple of 4
        if len % 4 != 0 {
            len += 4 - (len % 4);
        }

        // read string
        let mut bytes: Vec<u8> = vec![0; len as usize];
        self.reader.read_exact(&mut bytes)?;

        // convert to Rust string
        Ok(bytes2string(&bytes))
    }

    /// Read a string with one useful 8byte header from `XdrFile`.
    /// This is used for the body of the tpr files version >= 119.
    pub(super) fn read_string_8byte(&mut self) -> Result<String, Error> {
        // get length of the string
        let len = self.read_u64()?;

        // read string
        let mut bytes: Vec<u8> = vec![0; len as usize];
        self.reader.read_exact(&mut bytes)?;

        // convert to Rust string
        Ok(bytes2string(&bytes))
    }

    /// Read a string from the body of the tpr file.
    /// This calls either `read_string_4byte` or `read_string_8byte` depending on the
    /// version of the tpr file.
    #[inline(always)]
    pub(super) fn read_string_body(&mut self, tpr_version: i32) -> Result<String, Error> {
        if tpr_version < 119 {
            self.read_string_4byte()
        } else {
            self.read_string_8byte()
        }
    }

    /// Read a short unsigned integer from the body of the tpr file.
    /// This calls either `read_u32` or `read_u16` depending on the version of the tpr file.
    #[inline(always)]
    pub(super) fn read_ushort_body(&mut self, tpr_version: i32) -> Result<u32, Error> {
        if tpr_version < 119 {
            self.read_u32()
        } else {
            Ok(self.read_u16()? as u32)
        }
    }

    /// Read an unsigned char from the body of the tpr file.
    /// This calls either `read_u32` or `read_u8` depending on the version of the tpr file.
    #[inline(always)]
    pub(super) fn read_uchar_body(&mut self, tpr_version: i32) -> Result<u32, Error> {
        if tpr_version < 119 {
            self.read_u32()
        } else {
            Ok(self.read_u8()? as u32)
        }
    }

    /// Read a boolean value from the body of the tpr file.
    /// This calls either `read_u32` or `read_u8` depending on the version of the tpr file.
    #[inline(always)]
    pub(super) fn read_bool_body(&mut self, tpr_version: i32) -> Result<bool, Error> {
        Ok(self.read_uchar_body(tpr_version)? != 0)
    }

    /// Skip multiple unsigned char (or boolean) values from the body of the tpr file.
    /// This skips either `4n` or `n` bytes, depending on the version of the tpr file.
    #[inline(always)]
    pub(super) fn skip_multiple_uchars_body(
        &mut self,
        tpr_version: i32,
        n_uchars: i64,
    ) -> Result<(), Error> {
        if tpr_version < 119 {
            self.jump(4 * n_uchars)
        } else {
            self.jump(n_uchars)
        }
    }
}

/// Convert bytes to Rust string.
///
/// ## Panic
/// Panics if the conversion is not possible.
#[inline(always)]
fn bytes2string(bytes: &[u8]) -> String {
    let end = bytes.iter().position(|&b| b == 0).unwrap_or(bytes.len());

    std::str::from_utf8(&bytes[..end])
        .expect("FATAL MINITPR ERROR | bytes2string | Could not construct String from bytes.")
        .to_owned()
}
