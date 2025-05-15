// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024-2025 Ladislav Bartos

mod expected_values;

#[macro_use]
mod test_utilities {
    use float_cmp::assert_approx_eq;
    use minitpr::Atom;

    macro_rules! atom {
        ($atom_name:expr, $atom_number:expr, $residue_name:expr, $residue_number:expr, $mass:expr, $charge:expr, $element:expr, $position:expr, $velocity:expr, $force:expr) => {
            Atom {
                atom_name: $atom_name.to_owned(),
                atom_number: $atom_number,
                residue_name: $residue_name.to_owned(),
                residue_number: $residue_number,
                mass: $mass,
                charge: $charge,
                element: $element,
                position: $position,
                velocity: $velocity,
                force: $force,
            }
        };
    }

    macro_rules! bond {
        ($atom1:expr, $atom2:expr) => {
            Bond {
                atom1: $atom1,
                atom2: $atom2,
            }
        };
    }

    fn test_eq_coordinate(c1: &Option<[f64; 3]>, c2: &Option<[f64; 3]>) {
        match (c1, c2) {
            (None, None) => (),
            (Some(_), None) | (None, Some(_)) => {
                panic!("Coordinates do not match: {:?} != {:?}", c1, c2)
            }
            (Some(x), Some(y)) => {
                assert_approx_eq!(f64, x[0], y[0], epsilon = 0.001);
                assert_approx_eq!(f64, x[1], y[1], epsilon = 0.001);
                assert_approx_eq!(f64, x[2], y[2], epsilon = 0.001);
            }
        }
    }

    pub(super) fn test_eq_atom(atom: &Atom, expected: &Atom) {
        assert_eq!(atom.atom_name, expected.atom_name);
        assert_eq!(atom.atom_number, expected.atom_number);
        assert_eq!(atom.residue_name, expected.residue_name);
        assert_eq!(atom.residue_number, expected.residue_number);
        assert_approx_eq!(f64, atom.mass, expected.mass, epsilon = 0.000001);
        assert_approx_eq!(f64, atom.charge, expected.charge, epsilon = 0.000001);
        assert_eq!(atom.element, expected.element);
        test_eq_coordinate(&atom.position, &expected.position);
        test_eq_coordinate(&atom.velocity, &expected.velocity);
        test_eq_coordinate(&atom.force, &expected.force);
    }
}

#[cfg(test)]
mod tests {
    use super::test_utilities::*;
    use minitpr::{Atom, Bond, Element, Precision, TprFile};

    use float_cmp::assert_approx_eq;

    fn test_eq_small_cg(tpr: &TprFile, intermolecular: bool) {
        let header = &tpr.header;

        assert_eq!(header.precision, Precision::Single);
        assert_eq!(header.file_tag, "release");
        assert_eq!(header.n_atoms, 77);
        assert_eq!(header.n_coupling_groups, 1);
        assert_eq!(header.fep_state, 0);
        assert_approx_eq!(f64, header.lambda, 0.0);
        assert!(header.has_input_record);
        assert!(header.has_topology);
        assert!(header.has_positions);
        assert!(header.has_velocities);
        assert!(!header.has_forces);
        assert!(header.has_box);

        assert_eq!(&tpr.system_name, "Membrane");

        let simbox = tpr.simbox.as_ref().unwrap();

        assert_approx_eq!(f64, simbox.simbox[0][0], 9.2122, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[0][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[0][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox[1][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[1][1], 9.2122, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[1][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox[2][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[2][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[2][2], 11.3344, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_rel[0][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[0][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[0][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_rel[1][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[1][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[1][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_rel[2][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[2][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[2][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_v[0][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[0][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[0][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_v[1][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[1][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[1][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_v[2][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[2][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[2][2], 0.0000, epsilon = 0.000001);

        let topology = &tpr.topology;
        let atoms = &topology.atoms;

        assert_eq!(atoms.len(), 77);

        // expected positions of atoms
        let epos = crate::expected_values::CG_EXPECTED_POSITIONS;

        // expected velocity
        let evel = Some([0.0, 0.0, 0.0]);

        let expected = vec![
            // Peptide
            atom!("BB", 1, "LEU", 1, 72.0, 1.0, None, epos[0], evel, None),
            atom!("SC1", 2, "LEU", 1, 54.0, 0.0, None, epos[1], evel, None),
            atom!("BB", 3, "SER", 2, 72.0, 0.0, None, epos[2], evel, None),
            atom!("SC1", 4, "SER", 2, 36.0, 0.0, None, epos[3], evel, None),
            atom!("BB", 5, "SER", 3, 72.0, 0.0, None, epos[4], evel, None),
            atom!("SC1", 6, "SER", 3, 36.0, 0.0, None, epos[5], evel, None),
            atom!("BB", 7, "LEU", 4, 72.0, 0.0, None, epos[6], evel, None),
            atom!("SC1", 8, "LEU", 4, 54.0, 0.0, None, epos[7], evel, None),
            atom!("BB", 9, "LEU", 5, 72.0, 0.0, None, epos[8], evel, None),
            atom!("SC1", 10, "LEU", 5, 54.0, 0.0, None, epos[9], evel, None),
            atom!("BB", 11, "SER", 6, 72.0, 0.0, None, epos[10], evel, None),
            atom!("SC1", 12, "SER", 6, 36.0, 0.0, None, epos[11], evel, None),
            atom!("BB", 13, "LEU", 7, 72.0, 0.0, None, epos[12], evel, None),
            atom!("SC1", 14, "LEU", 7, 54.0, 0.0, None, epos[13], evel, None),
            atom!("BB", 15, "LEU", 8, 72.0, 0.0, None, epos[14], evel, None),
            atom!("SC1", 16, "LEU", 8, 54.0, 0.0, None, epos[15], evel, None),
            atom!("BB", 17, "SER", 9, 72.0, 0.0, None, epos[16], evel, None),
            atom!("SC1", 18, "SER", 9, 36.0, 0.0, None, epos[17], evel, None),
            atom!("BB", 19, "SER", 10, 72.0, 0.0, None, epos[18], evel, None),
            atom!("SC1", 20, "SER", 10, 36.0, 0.0, None, epos[19], evel, None),
            atom!("BB", 21, "LEU", 11, 72.0, 0.0, None, epos[20], evel, None),
            atom!("SC1", 22, "LEU", 11, 54.0, 0.0, None, epos[21], evel, None),
            atom!("BB", 23, "LEU", 12, 72.0, 0.0, None, epos[22], evel, None),
            atom!("SC1", 24, "LEU", 12, 54.0, 0.0, None, epos[23], evel, None),
            atom!("BB", 25, "SER", 13, 72.0, 0.0, None, epos[24], evel, None),
            atom!("SC1", 26, "SER", 13, 36.0, 0.0, None, epos[25], evel, None),
            atom!("BB", 27, "LEU", 14, 72.0, 0.0, None, epos[26], evel, None),
            atom!("SC1", 28, "LEU", 14, 54.0, 0.0, None, epos[27], evel, None),
            atom!("BB", 29, "LEU", 15, 72.0, 0.0, None, epos[28], evel, None),
            atom!("SC1", 30, "LEU", 15, 54.0, 0.0, None, epos[29], evel, None),
            atom!("BB", 31, "SER", 16, 72.0, 0.0, None, epos[30], evel, None),
            atom!("SC1", 32, "SER", 16, 36.0, 0.0, None, epos[31], evel, None),
            atom!("BB", 33, "SER", 17, 72.0, 0.0, None, epos[32], evel, None),
            atom!("SC1", 34, "SER", 17, 36.0, 0.0, None, epos[33], evel, None),
            atom!("BB", 35, "LEU", 18, 72.0, 0.0, None, epos[34], evel, None),
            atom!("SC1", 36, "LEU", 18, 54.0, 0.0, None, epos[35], evel, None),
            atom!("BB", 37, "LEU", 19, 72.0, 0.0, None, epos[36], evel, None),
            atom!("SC1", 38, "LEU", 19, 54.0, 0.0, None, epos[37], evel, None),
            atom!("BB", 39, "SER", 20, 72.0, 0.0, None, epos[38], evel, None),
            atom!("SC1", 40, "SER", 20, 36.0, 0.0, None, epos[39], evel, None),
            atom!("BB", 41, "LEU", 21, 72.0, 0.0, None, epos[40], evel, None),
            atom!("SC1", 42, "LEU", 21, 54.0, 0.0, None, epos[41], evel, None),
            // POPC #1
            atom!("NC3", 43, "POPC", 22, 72.0, 1.0, None, epos[42], evel, None),
            atom!("PO4", 44, "POPC", 22, 72.0, -1.0, None, epos[43], evel, None),
            atom!("GL1", 45, "POPC", 22, 54.0, 0.0, None, epos[44], evel, None),
            atom!("GL2", 46, "POPC", 22, 72.0, 0.0, None, epos[45], evel, None),
            atom!("C1A", 47, "POPC", 22, 72.0, 0.0, None, epos[46], evel, None),
            atom!("D2A", 48, "POPC", 22, 72.0, 0.0, None, epos[47], evel, None),
            atom!("C3A", 49, "POPC", 22, 72.0, 0.0, None, epos[48], evel, None),
            atom!("C4A", 50, "POPC", 22, 72.0, 0.0, None, epos[49], evel, None),
            atom!("C1B", 51, "POPC", 22, 72.0, 0.0, None, epos[50], evel, None),
            atom!("C2B", 52, "POPC", 22, 72.0, 0.0, None, epos[51], evel, None),
            atom!("C3B", 53, "POPC", 22, 72.0, 0.0, None, epos[52], evel, None),
            atom!("C4B", 54, "POPC", 22, 72.0, 0.0, None, epos[53], evel, None),
            // POPC #2
            atom!("NC3", 55, "POPC", 23, 72.0, 1.0, None, epos[54], evel, None),
            atom!("PO4", 56, "POPC", 23, 72.0, -1.0, None, epos[55], evel, None),
            atom!("GL1", 57, "POPC", 23, 54.0, 0.0, None, epos[56], evel, None),
            atom!("GL2", 58, "POPC", 23, 72.0, 0.0, None, epos[57], evel, None),
            atom!("C1A", 59, "POPC", 23, 72.0, 0.0, None, epos[58], evel, None),
            atom!("D2A", 60, "POPC", 23, 72.0, 0.0, None, epos[59], evel, None),
            atom!("C3A", 61, "POPC", 23, 72.0, 0.0, None, epos[60], evel, None),
            atom!("C4A", 62, "POPC", 23, 72.0, 0.0, None, epos[61], evel, None),
            atom!("C1B", 63, "POPC", 23, 72.0, 0.0, None, epos[62], evel, None),
            atom!("C2B", 64, "POPC", 23, 72.0, 0.0, None, epos[63], evel, None),
            atom!("C3B", 65, "POPC", 23, 72.0, 0.0, None, epos[64], evel, None),
            atom!("C4B", 66, "POPC", 23, 72.0, 0.0, None, epos[65], evel, None),
            // water (10x)
            atom!("W", 67, "W", 24, 72.0, 0.0, None, epos[66], evel, None),
            atom!("W", 68, "W", 25, 72.0, 0.0, None, epos[67], evel, None),
            atom!("W", 69, "W", 26, 72.0, 0.0, None, epos[68], evel, None),
            atom!("W", 70, "W", 27, 72.0, 0.0, None, epos[69], evel, None),
            atom!("W", 71, "W", 28, 72.0, 0.0, None, epos[70], evel, None),
            atom!("W", 72, "W", 29, 72.0, 0.0, None, epos[71], evel, None),
            atom!("W", 73, "W", 30, 72.0, 0.0, None, epos[72], evel, None),
            atom!("W", 74, "W", 31, 72.0, 0.0, None, epos[73], evel, None),
            atom!("W", 75, "W", 32, 72.0, 0.0, None, epos[74], evel, None),
            atom!("W", 76, "W", 33, 72.0, 0.0, None, epos[75], evel, None),
            // CL ion
            atom!("CL-", 77, "ION", 34, 35.453, -1.0, None, epos[76], evel, None),
        ];

        for (a, e) in atoms.iter().zip(expected.iter()) {
            test_eq_atom(a, e);
        }

        let bonds = &topology.bonds;

        let mut expected_bonds = vec![
            bond!(0, 1),
            bond!(2, 3),
            bond!(4, 5),
            bond!(6, 7),
            bond!(8, 9),
            bond!(10, 11),
            bond!(12, 13),
            bond!(14, 15),
            bond!(16, 17),
            bond!(18, 19),
            bond!(20, 21),
            bond!(22, 23),
            bond!(24, 25),
            bond!(26, 27),
            bond!(28, 29),
            bond!(30, 31),
            bond!(32, 33),
            bond!(34, 35),
            bond!(36, 37),
            bond!(38, 39),
            bond!(40, 41),
            bond!(42, 43),
            bond!(43, 44),
            bond!(44, 45),
            bond!(44, 46),
            bond!(46, 47),
            bond!(47, 48),
            bond!(48, 49),
            bond!(45, 50),
            bond!(50, 51),
            bond!(51, 52),
            bond!(52, 53),
            bond!(54, 55),
            bond!(55, 56),
            bond!(56, 57),
            bond!(56, 58),
            bond!(58, 59),
            bond!(59, 60),
            bond!(60, 61),
            bond!(57, 62),
            bond!(62, 63),
            bond!(63, 64),
            bond!(64, 65),
            bond!(0, 2),
            bond!(2, 4),
            bond!(4, 6),
            bond!(6, 8),
            bond!(8, 10),
            bond!(10, 12),
            bond!(12, 14),
            bond!(14, 16),
            bond!(16, 18),
            bond!(18, 20),
            bond!(20, 22),
            bond!(22, 24),
            bond!(24, 26),
            bond!(26, 28),
            bond!(28, 30),
            bond!(30, 32),
            bond!(32, 34),
            bond!(34, 36),
            bond!(36, 38),
            bond!(38, 40),
        ];

        if intermolecular {
            expected_bonds.extend(vec![bond!(28, 52), bond!(64, 67), bond!(73, 74)]);
        }

        assert_eq!(bonds.len(), expected_bonds.len());

        for e in expected_bonds.iter() {
            assert!(bonds.contains(e));
        }
    }

    #[test]
    fn small_cg_5() {
        let tpr = TprFile::parse("tests/test_files/small_cg_5.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 5.1.4"));
        assert_eq!(header.tpr_version, 103);
        assert_eq!(header.tpr_generation, 26);
        assert!(header.body_size.is_none());

        test_eq_small_cg(&tpr, false);
    }

    #[test]
    fn small_cg_2016() {
        let tpr = TprFile::parse("tests/test_files/small_cg_2016.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2016.4"));
        assert_eq!(header.tpr_version, 110);
        assert_eq!(header.tpr_generation, 26);
        assert!(header.body_size.is_none());

        test_eq_small_cg(&tpr, false);
    }

    #[test]
    fn small_cg_2021() {
        let tpr = TprFile::parse("tests/test_files/small_cg_2021.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2021.4"));
        assert_eq!(header.tpr_version, 122);
        assert_eq!(header.tpr_generation, 28);
        assert_eq!(header.body_size.unwrap(), 24909);

        test_eq_small_cg(&tpr, false);
    }

    #[test]
    fn small_cg_5_intermolecular() {
        let tpr = TprFile::parse("tests/test_files/small_cg_5_intermolecular.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 5.1.4"));
        assert_eq!(header.tpr_version, 103);
        assert_eq!(header.tpr_generation, 26);
        assert!(header.body_size.is_none());

        test_eq_small_cg(&tpr, true);
    }

    #[test]
    fn small_cg_2016_intermolecular() {
        let tpr = TprFile::parse("tests/test_files/small_cg_2016_intermolecular.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2016.4"));
        assert_eq!(header.tpr_version, 110);
        assert_eq!(header.tpr_generation, 26);
        assert!(header.body_size.is_none());

        test_eq_small_cg(&tpr, true);
    }

    #[test]
    fn small_cg_2021_intermolecular() {
        let tpr = TprFile::parse("tests/test_files/small_cg_2021_intermolecular.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2021.4"));
        assert_eq!(header.tpr_version, 122);
        assert_eq!(header.tpr_generation, 28);
        assert_eq!(header.body_size.unwrap(), 25429);

        test_eq_small_cg(&tpr, true);
    }

    #[test]
    fn empty_fail() {
        assert!(TprFile::parse("tests/test_files/empty.tpr").is_err());
    }

    enum GmxVersion {
        Gromacs5,
        Gromacs2016,
        Gromacs2021,
    }

    fn test_eq_small_aa(tpr: &TprFile, intermolecular: bool, version: GmxVersion) {
        let header = &tpr.header;

        assert_eq!(header.precision, Precision::Single);
        assert_eq!(header.file_tag, "release");
        assert_eq!(header.n_atoms, 182);
        assert_eq!(header.n_coupling_groups, 3);
        assert_eq!(header.fep_state, 0);
        assert_approx_eq!(f64, header.lambda, 0.0);
        assert!(header.has_input_record);
        assert!(header.has_topology);
        assert!(header.has_positions);
        assert!(header.has_velocities);
        assert!(!header.has_forces);
        assert!(header.has_box);

        assert_eq!(&tpr.system_name, "Protein");

        let simbox = tpr.simbox.as_ref().unwrap();

        assert_approx_eq!(f64, simbox.simbox[0][0], 10.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[0][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[0][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox[1][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[1][1], 10.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[1][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox[2][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[2][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[2][2], 10.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_rel[0][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[0][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[0][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_rel[1][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[1][1], 1.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[1][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_rel[2][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[2][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[2][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_v[0][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[0][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[0][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_v[1][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[1][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[1][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_v[2][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[2][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[2][2], 0.0000, epsilon = 0.000001);

        let topology = &tpr.topology;
        let atoms = &topology.atoms;

        assert_eq!(atoms.len(), 182);

        // expected positions
        let epos = crate::expected_values::AA_EXPECTED_POSITIONS;

        // expected velocities
        let evel = match (version, intermolecular) {
            (GmxVersion::Gromacs2021, false) => {
                crate::expected_values::AA_GMX2021_EXPECTED_VELOCITIES
            }
            (GmxVersion::Gromacs2021, true) => {
                crate::expected_values::AA_GMX2021_EXPECTED_VELOCITIES_INTERMOLECULAR
            }
            (GmxVersion::Gromacs2016, false) => {
                crate::expected_values::AA_GMX2016_EXPECTED_VELOCITIES
            }
            (GmxVersion::Gromacs2016, true) => {
                crate::expected_values::AA_GMX2016_EXPECTED_VELOCITIES_INTERMOLECULAR
            }
            (GmxVersion::Gromacs5, false) => crate::expected_values::AA_GMX5_EXPECTED_VELOCITIES,
            (GmxVersion::Gromacs5, true) => {
                crate::expected_values::AA_GMX5_EXPECTED_VELOCITIES_INTERMOLECULAR
            }
        };

        let expected = vec![
            // dipeptide
            atom!(
                "N",
                1,
                "LEU",
                1,
                14.01,
                0.101,
                Some(Element::N),
                epos[0],
                evel[0],
                None
            ),
            atom!(
                "H1",
                2,
                "LEU",
                1,
                1.008,
                0.2148,
                Some(Element::H),
                epos[1],
                evel[1],
                None
            ),
            atom!(
                "H2",
                3,
                "LEU",
                1,
                1.008,
                0.2148,
                Some(Element::H),
                epos[2],
                evel[2],
                None
            ),
            atom!(
                "H3",
                4,
                "LEU",
                1,
                1.008,
                0.2148,
                Some(Element::H),
                epos[3],
                evel[3],
                None
            ),
            atom!(
                "CA",
                5,
                "LEU",
                1,
                12.01,
                0.0104,
                Some(Element::C),
                epos[4],
                evel[4],
                None
            ),
            atom!(
                "HA",
                6,
                "LEU",
                1,
                1.008,
                0.1053,
                Some(Element::H),
                epos[5],
                evel[5],
                None
            ),
            atom!(
                "CB",
                7,
                "LEU",
                1,
                12.01,
                -0.0244,
                Some(Element::C),
                epos[6],
                evel[6],
                None
            ),
            atom!(
                "HB1",
                8,
                "LEU",
                1,
                1.008,
                0.0256,
                Some(Element::H),
                epos[7],
                evel[7],
                None
            ),
            atom!(
                "HB2",
                9,
                "LEU",
                1,
                1.008,
                0.0256,
                Some(Element::H),
                epos[8],
                evel[8],
                None
            ),
            atom!(
                "CG",
                10,
                "LEU",
                1,
                12.01,
                0.3421,
                Some(Element::C),
                epos[9],
                evel[9],
                None
            ),
            atom!(
                "HG",
                11,
                "LEU",
                1,
                1.008,
                -0.038,
                Some(Element::H),
                epos[10],
                evel[10],
                None
            ),
            atom!(
                "CD1",
                12,
                "LEU",
                1,
                12.01,
                -0.4106,
                Some(Element::C),
                epos[11],
                evel[11],
                None
            ),
            atom!(
                "HD11",
                13,
                "LEU",
                1,
                1.008,
                0.098,
                Some(Element::H),
                epos[12],
                evel[12],
                None
            ),
            atom!(
                "HD12",
                14,
                "LEU",
                1,
                1.008,
                0.098,
                Some(Element::H),
                epos[13],
                evel[13],
                None
            ),
            atom!(
                "HD13",
                15,
                "LEU",
                1,
                1.008,
                0.098,
                Some(Element::H),
                epos[14],
                evel[14],
                None
            ),
            atom!(
                "CD2",
                16,
                "LEU",
                1,
                12.01,
                -0.4104,
                Some(Element::C),
                epos[15],
                evel[15],
                None
            ),
            atom!(
                "HD21",
                17,
                "LEU",
                1,
                1.008,
                0.098,
                Some(Element::H),
                epos[16],
                evel[16],
                None
            ),
            atom!(
                "HD22",
                18,
                "LEU",
                1,
                1.008,
                0.098,
                Some(Element::H),
                epos[17],
                evel[17],
                None
            ),
            atom!(
                "HD23",
                19,
                "LEU",
                1,
                1.008,
                0.098,
                Some(Element::H),
                epos[18],
                evel[18],
                None
            ),
            atom!(
                "C",
                20,
                "LEU",
                1,
                12.01,
                0.6123,
                Some(Element::C),
                epos[19],
                evel[19],
                None
            ),
            atom!(
                "O",
                21,
                "LEU",
                1,
                16.00,
                -0.5713,
                Some(Element::O),
                epos[20],
                evel[20],
                None
            ),
            atom!(
                "N",
                22,
                "LYS",
                2,
                14.01,
                -0.3481,
                Some(Element::N),
                epos[21],
                evel[21],
                None
            ),
            atom!(
                "H",
                23,
                "LYS",
                2,
                1.008,
                0.2764,
                Some(Element::H),
                epos[22],
                evel[22],
                None
            ),
            atom!(
                "CA",
                24,
                "LYS",
                2,
                12.01,
                -0.2903,
                Some(Element::C),
                epos[23],
                evel[23],
                None
            ),
            atom!(
                "HA",
                25,
                "LYS",
                2,
                1.008,
                0.1438,
                Some(Element::H),
                epos[24],
                evel[24],
                None
            ),
            atom!(
                "CB",
                26,
                "LYS",
                2,
                12.01,
                -0.0538,
                Some(Element::C),
                epos[25],
                evel[25],
                None
            ),
            atom!(
                "HB1",
                27,
                "LYS",
                2,
                1.008,
                0.0482,
                Some(Element::H),
                epos[26],
                evel[26],
                None
            ),
            atom!(
                "HB2",
                28,
                "LYS",
                2,
                1.008,
                0.0482,
                Some(Element::H),
                epos[27],
                evel[27],
                None
            ),
            atom!(
                "CG",
                29,
                "LYS",
                2,
                12.01,
                0.0227,
                Some(Element::C),
                epos[28],
                evel[28],
                None
            ),
            atom!(
                "HG1",
                30,
                "LYS",
                2,
                1.008,
                0.0134,
                Some(Element::H),
                epos[29],
                evel[29],
                None
            ),
            atom!(
                "HG2",
                31,
                "LYS",
                2,
                1.008,
                0.0134,
                Some(Element::H),
                epos[30],
                evel[30],
                None
            ),
            atom!(
                "CD",
                32,
                "LYS",
                2,
                12.01,
                -0.0392,
                Some(Element::C),
                epos[31],
                evel[31],
                None
            ),
            atom!(
                "HD1",
                33,
                "LYS",
                2,
                1.008,
                0.0611,
                Some(Element::H),
                epos[32],
                evel[32],
                None
            ),
            atom!(
                "HD2",
                34,
                "LYS",
                2,
                1.008,
                0.0611,
                Some(Element::H),
                epos[33],
                evel[33],
                None
            ),
            atom!(
                "CE",
                35,
                "LYS",
                2,
                12.01,
                -0.0176,
                Some(Element::C),
                epos[34],
                evel[34],
                None
            ),
            atom!(
                "HE1",
                36,
                "LYS",
                2,
                1.008,
                0.1121,
                Some(Element::H),
                epos[35],
                evel[35],
                None
            ),
            atom!(
                "HE2",
                37,
                "LYS",
                2,
                1.008,
                0.1121,
                Some(Element::H),
                epos[36],
                evel[36],
                None
            ),
            atom!(
                "NZ",
                38,
                "LYS",
                2,
                14.01,
                -0.3741,
                Some(Element::N),
                epos[37],
                evel[37],
                None
            ),
            atom!(
                "HZ1",
                39,
                "LYS",
                2,
                1.008,
                0.3374,
                Some(Element::H),
                epos[38],
                evel[38],
                None
            ),
            atom!(
                "HZ2",
                40,
                "LYS",
                2,
                1.008,
                0.3374,
                Some(Element::H),
                epos[39],
                evel[39],
                None
            ),
            atom!(
                "HZ3",
                41,
                "LYS",
                2,
                1.008,
                0.3374,
                Some(Element::H),
                epos[40],
                evel[40],
                None
            ),
            atom!(
                "C",
                42,
                "LYS",
                2,
                12.01,
                0.8488,
                Some(Element::C),
                epos[41],
                evel[41],
                None
            ),
            atom!(
                "OC1",
                43,
                "LYS",
                2,
                16.00,
                -0.8252,
                Some(Element::O),
                epos[42],
                evel[42],
                None
            ),
            atom!(
                "OC2",
                44,
                "LYS",
                2,
                16.00,
                -0.8252,
                Some(Element::O),
                epos[43],
                evel[43],
                None
            ),
            // POPC
            atom!(
                "N",
                45,
                "POPC",
                3,
                14.007,
                0.20,
                Some(Element::N),
                epos[44],
                evel[44],
                None
            ),
            atom!(
                "C12",
                46,
                "POPC",
                3,
                12.011,
                -0.20,
                Some(Element::C),
                epos[45],
                evel[45],
                None
            ),
            atom!(
                "C13",
                47,
                "POPC",
                3,
                12.011,
                -0.38,
                Some(Element::C),
                epos[46],
                evel[46],
                None
            ),
            atom!(
                "C14",
                48,
                "POPC",
                3,
                12.011,
                -0.38,
                Some(Element::C),
                epos[47],
                evel[47],
                None
            ),
            atom!(
                "C15",
                49,
                "POPC",
                3,
                12.011,
                -0.38,
                Some(Element::C),
                epos[48],
                evel[48],
                None
            ),
            atom!(
                "H12A",
                50,
                "POPC",
                3,
                1.008,
                0.09,
                Some(Element::H),
                epos[49],
                evel[49],
                None
            ),
            atom!(
                "H12B",
                51,
                "POPC",
                3,
                1.008,
                0.09,
                Some(Element::H),
                epos[50],
                evel[50],
                None
            ),
            atom!(
                "H13A",
                52,
                "POPC",
                3,
                1.008,
                0.19,
                Some(Element::H),
                epos[51],
                evel[51],
                None
            ),
            atom!(
                "H13B",
                53,
                "POPC",
                3,
                1.008,
                0.19,
                Some(Element::H),
                epos[52],
                evel[52],
                None
            ),
            atom!(
                "H13C",
                54,
                "POPC",
                3,
                1.008,
                0.19,
                Some(Element::H),
                epos[53],
                evel[53],
                None
            ),
            atom!(
                "H14A",
                55,
                "POPC",
                3,
                1.008,
                0.19,
                Some(Element::H),
                epos[54],
                evel[54],
                None
            ),
            atom!(
                "H14B",
                56,
                "POPC",
                3,
                1.008,
                0.19,
                Some(Element::H),
                epos[55],
                evel[55],
                None
            ),
            atom!(
                "H14C",
                57,
                "POPC",
                3,
                1.008,
                0.19,
                Some(Element::H),
                epos[56],
                evel[56],
                None
            ),
            atom!(
                "H15A",
                58,
                "POPC",
                3,
                1.008,
                0.19,
                Some(Element::H),
                epos[57],
                evel[57],
                None
            ),
            atom!(
                "H15B",
                59,
                "POPC",
                3,
                1.008,
                0.19,
                Some(Element::H),
                epos[58],
                evel[58],
                None
            ),
            atom!(
                "H15C",
                60,
                "POPC",
                3,
                1.008,
                0.19,
                Some(Element::H),
                epos[59],
                evel[59],
                None
            ),
            atom!(
                "C11",
                61,
                "POPC",
                3,
                12.011,
                0.17,
                Some(Element::C),
                epos[60],
                evel[60],
                None
            ),
            atom!(
                "H11A",
                62,
                "POPC",
                3,
                1.008,
                0.03,
                Some(Element::H),
                epos[61],
                evel[61],
                None
            ),
            atom!(
                "H11B",
                63,
                "POPC",
                3,
                1.008,
                0.03,
                Some(Element::H),
                epos[62],
                evel[62],
                None
            ),
            atom!(
                "P",
                64,
                "POPC",
                3,
                30.974,
                1.58,
                Some(Element::P),
                epos[63],
                evel[63],
                None
            ),
            atom!(
                "O13",
                65,
                "POPC",
                3,
                15.9994,
                -0.86,
                Some(Element::O),
                epos[64],
                evel[64],
                None
            ),
            atom!(
                "O14",
                66,
                "POPC",
                3,
                15.9994,
                -0.86,
                Some(Element::O),
                epos[65],
                evel[65],
                None
            ),
            atom!(
                "O12",
                67,
                "POPC",
                3,
                15.9994,
                -0.49,
                Some(Element::O),
                epos[66],
                evel[66],
                None
            ),
            atom!(
                "O11",
                68,
                "POPC",
                3,
                15.9994,
                -0.49,
                Some(Element::O),
                epos[67],
                evel[67],
                None
            ),
            atom!(
                "C1",
                69,
                "POPC",
                3,
                12.011,
                -0.11,
                Some(Element::C),
                epos[68],
                evel[68],
                None
            ),
            atom!(
                "HA",
                70,
                "POPC",
                3,
                1.008,
                0.07,
                Some(Element::H),
                epos[69],
                evel[69],
                None
            ),
            atom!(
                "HB",
                71,
                "POPC",
                3,
                1.008,
                0.07,
                Some(Element::H),
                epos[70],
                evel[70],
                None
            ),
            atom!(
                "C2",
                72,
                "POPC",
                3,
                12.011,
                0.48,
                Some(Element::C),
                epos[71],
                evel[71],
                None
            ),
            atom!(
                "HS",
                73,
                "POPC",
                3,
                1.008,
                0.04,
                Some(Element::H),
                epos[72],
                evel[72],
                None
            ),
            atom!(
                "O21",
                74,
                "POPC",
                3,
                15.9994,
                -0.47,
                Some(Element::O),
                epos[73],
                evel[73],
                None
            ),
            atom!(
                "C21",
                75,
                "POPC",
                3,
                12.011,
                0.79,
                Some(Element::C),
                epos[74],
                evel[74],
                None
            ),
            atom!(
                "O22",
                76,
                "POPC",
                3,
                15.9994,
                -0.65,
                Some(Element::O),
                epos[75],
                evel[75],
                None
            ),
            atom!(
                "C22",
                77,
                "POPC",
                3,
                12.011,
                -0.06,
                Some(Element::C),
                epos[76],
                evel[76],
                None
            ),
            atom!(
                "H2R",
                78,
                "POPC",
                3,
                1.008,
                0.03,
                Some(Element::H),
                epos[77],
                evel[77],
                None
            ),
            atom!(
                "H2S",
                79,
                "POPC",
                3,
                1.008,
                0.03,
                Some(Element::H),
                epos[78],
                evel[78],
                None
            ),
            atom!(
                "C3",
                80,
                "POPC",
                3,
                12.011,
                0.13,
                Some(Element::C),
                epos[79],
                evel[79],
                None
            ),
            atom!(
                "HX",
                81,
                "POPC",
                3,
                1.008,
                0.06,
                Some(Element::H),
                epos[80],
                evel[80],
                None
            ),
            atom!(
                "HY",
                82,
                "POPC",
                3,
                1.008,
                0.06,
                Some(Element::H),
                epos[81],
                evel[81],
                None
            ),
            atom!(
                "O31",
                83,
                "POPC",
                3,
                15.9994,
                -0.47,
                Some(Element::O),
                epos[82],
                evel[82],
                None
            ),
            atom!(
                "C31",
                84,
                "POPC",
                3,
                12.011,
                0.79,
                Some(Element::C),
                epos[83],
                evel[83],
                None
            ),
            atom!(
                "O32",
                85,
                "POPC",
                3,
                15.9994,
                -0.65,
                Some(Element::O),
                epos[84],
                evel[84],
                None
            ),
            atom!(
                "C32",
                86,
                "POPC",
                3,
                12.011,
                -0.06,
                Some(Element::C),
                epos[85],
                evel[85],
                None
            ),
            atom!(
                "H2X",
                87,
                "POPC",
                3,
                1.008,
                0.03,
                Some(Element::H),
                epos[86],
                evel[86],
                None
            ),
            atom!(
                "H2Y",
                88,
                "POPC",
                3,
                1.008,
                0.03,
                Some(Element::H),
                epos[87],
                evel[87],
                None
            ),
            atom!(
                "C23",
                89,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[88],
                evel[88],
                None
            ),
            atom!(
                "H3R",
                90,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[89],
                evel[89],
                None
            ),
            atom!(
                "H3S",
                91,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[90],
                evel[90],
                None
            ),
            atom!(
                "C24",
                92,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[91],
                evel[91],
                None
            ),
            atom!(
                "H4R",
                93,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[92],
                evel[92],
                None
            ),
            atom!(
                "H4S",
                94,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[93],
                evel[93],
                None
            ),
            atom!(
                "C25",
                95,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[94],
                evel[94],
                None
            ),
            atom!(
                "H5R",
                96,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[95],
                evel[95],
                None
            ),
            atom!(
                "H5S",
                97,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[96],
                evel[96],
                None
            ),
            atom!(
                "C26",
                98,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[97],
                evel[97],
                None
            ),
            atom!(
                "H6R",
                99,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[98],
                evel[98],
                None
            ),
            atom!(
                "H6S",
                100,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[99],
                evel[99],
                None
            ),
            atom!(
                "C27",
                101,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[100],
                evel[100],
                None
            ),
            atom!(
                "H7R",
                102,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[101],
                evel[101],
                None
            ),
            atom!(
                "H7S",
                103,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[102],
                evel[102],
                None
            ),
            atom!(
                "C28",
                104,
                "POPC",
                3,
                12.011,
                0.03,
                Some(Element::C),
                epos[103],
                evel[103],
                None
            ),
            atom!(
                "H8R",
                105,
                "POPC",
                3,
                1.008,
                0.03,
                Some(Element::H),
                epos[104],
                evel[104],
                None
            ),
            atom!(
                "H8S",
                106,
                "POPC",
                3,
                1.008,
                0.03,
                Some(Element::H),
                epos[105],
                evel[105],
                None
            ),
            atom!(
                "C29",
                107,
                "POPC",
                3,
                12.011,
                -0.20,
                Some(Element::C),
                epos[106],
                evel[106],
                None
            ),
            atom!(
                "H91",
                108,
                "POPC",
                3,
                1.008,
                0.11,
                Some(Element::H),
                epos[107],
                evel[107],
                None
            ),
            atom!(
                "C210",
                109,
                "POPC",
                3,
                12.011,
                -0.20,
                Some(Element::C),
                epos[108],
                evel[108],
                None
            ),
            atom!(
                "H101",
                110,
                "POPC",
                3,
                1.008,
                0.11,
                Some(Element::H),
                epos[109],
                evel[109],
                None
            ),
            atom!(
                "C211",
                111,
                "POPC",
                3,
                12.011,
                0.03,
                Some(Element::C),
                epos[110],
                evel[110],
                None
            ),
            atom!(
                "H11R",
                112,
                "POPC",
                3,
                1.008,
                0.03,
                Some(Element::H),
                epos[111],
                evel[111],
                None
            ),
            atom!(
                "H11S",
                113,
                "POPC",
                3,
                1.008,
                0.03,
                Some(Element::H),
                epos[112],
                evel[112],
                None
            ),
            atom!(
                "C212",
                114,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[113],
                evel[113],
                None
            ),
            atom!(
                "H12R",
                115,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[114],
                evel[114],
                None
            ),
            atom!(
                "H12S",
                116,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[115],
                evel[115],
                None
            ),
            atom!(
                "C213",
                117,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[116],
                evel[116],
                None
            ),
            atom!(
                "H13R",
                118,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[117],
                evel[117],
                None
            ),
            atom!(
                "H13S",
                119,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[118],
                evel[118],
                None
            ),
            atom!(
                "C214",
                120,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[119],
                evel[119],
                None
            ),
            atom!(
                "H14R",
                121,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[120],
                evel[120],
                None
            ),
            atom!(
                "H14S",
                122,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[121],
                evel[121],
                None
            ),
            atom!(
                "C215",
                123,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[122],
                evel[122],
                None
            ),
            atom!(
                "H15R",
                124,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[123],
                evel[123],
                None
            ),
            atom!(
                "H15S",
                125,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[124],
                evel[124],
                None
            ),
            atom!(
                "C216",
                126,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[125],
                evel[125],
                None
            ),
            atom!(
                "H16R",
                127,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[126],
                evel[126],
                None
            ),
            atom!(
                "H16S",
                128,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[127],
                evel[127],
                None
            ),
            atom!(
                "C217",
                129,
                "POPC",
                3,
                12.011,
                0.047,
                Some(Element::C),
                epos[128],
                evel[128],
                None
            ),
            atom!(
                "H17R",
                130,
                "POPC",
                3,
                1.008,
                -0.007,
                Some(Element::H),
                epos[129],
                evel[129],
                None
            ),
            atom!(
                "H17S",
                131,
                "POPC",
                3,
                1.008,
                -0.007,
                Some(Element::H),
                epos[130],
                evel[130],
                None
            ),
            atom!(
                "C218",
                132,
                "POPC",
                3,
                12.011,
                -0.081,
                Some(Element::C),
                epos[131],
                evel[131],
                None
            ),
            atom!(
                "H18R",
                133,
                "POPC",
                3,
                1.008,
                0.016,
                Some(Element::H),
                epos[132],
                evel[132],
                None
            ),
            atom!(
                "H18S",
                134,
                "POPC",
                3,
                1.008,
                0.016,
                Some(Element::H),
                epos[133],
                evel[133],
                None
            ),
            atom!(
                "H18T",
                135,
                "POPC",
                3,
                1.008,
                0.016,
                Some(Element::H),
                epos[134],
                evel[134],
                None
            ),
            atom!(
                "C33",
                136,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[135],
                evel[135],
                None
            ),
            atom!(
                "H3X",
                137,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[136],
                evel[136],
                None
            ),
            atom!(
                "H3Y",
                138,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[137],
                evel[137],
                None
            ),
            atom!(
                "C34",
                139,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[138],
                evel[138],
                None
            ),
            atom!(
                "H4X",
                140,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[139],
                evel[139],
                None
            ),
            atom!(
                "H4Y",
                141,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[140],
                evel[140],
                None
            ),
            atom!(
                "C35",
                142,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[141],
                evel[141],
                None
            ),
            atom!(
                "H5X",
                143,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[142],
                evel[142],
                None
            ),
            atom!(
                "H5Y",
                144,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[143],
                evel[143],
                None
            ),
            atom!(
                "C36",
                145,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[144],
                evel[144],
                None
            ),
            atom!(
                "H6X",
                146,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[145],
                evel[145],
                None
            ),
            atom!(
                "H6Y",
                147,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[146],
                evel[146],
                None
            ),
            atom!(
                "C37",
                148,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[147],
                evel[147],
                None
            ),
            atom!(
                "H7X",
                149,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[148],
                evel[148],
                None
            ),
            atom!(
                "H7Y",
                150,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[149],
                evel[149],
                None
            ),
            atom!(
                "C38",
                151,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[150],
                evel[150],
                None
            ),
            atom!(
                "H8X",
                152,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[151],
                evel[151],
                None
            ),
            atom!(
                "H8Y",
                153,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[152],
                evel[152],
                None
            ),
            atom!(
                "C39",
                154,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[153],
                evel[153],
                None
            ),
            atom!(
                "H9X",
                155,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[154],
                evel[154],
                None
            ),
            atom!(
                "H9Y",
                156,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[155],
                evel[155],
                None
            ),
            atom!(
                "C310",
                157,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[156],
                evel[156],
                None
            ),
            atom!(
                "H10X",
                158,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[157],
                evel[157],
                None
            ),
            atom!(
                "H10Y",
                159,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[158],
                evel[158],
                None
            ),
            atom!(
                "C311",
                160,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[159],
                evel[159],
                None
            ),
            atom!(
                "H11X",
                161,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[160],
                evel[160],
                None
            ),
            atom!(
                "H11Y",
                162,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[161],
                evel[161],
                None
            ),
            atom!(
                "C312",
                163,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[162],
                evel[162],
                None
            ),
            atom!(
                "H12X",
                164,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[163],
                evel[163],
                None
            ),
            atom!(
                "H12Y",
                165,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[164],
                evel[164],
                None
            ),
            atom!(
                "C313",
                166,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[165],
                evel[165],
                None
            ),
            atom!(
                "H13X",
                167,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[166],
                evel[166],
                None
            ),
            atom!(
                "H13Y",
                168,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[167],
                evel[167],
                None
            ),
            atom!(
                "C314",
                169,
                "POPC",
                3,
                12.011,
                0.00,
                Some(Element::C),
                epos[168],
                evel[168],
                None
            ),
            atom!(
                "H14X",
                170,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[169],
                evel[169],
                None
            ),
            atom!(
                "H14Y",
                171,
                "POPC",
                3,
                1.008,
                0.00,
                Some(Element::H),
                epos[170],
                evel[170],
                None
            ),
            atom!(
                "C315",
                172,
                "POPC",
                3,
                12.011,
                0.047,
                Some(Element::C),
                epos[171],
                evel[171],
                None
            ),
            atom!(
                "H15X",
                173,
                "POPC",
                3,
                1.008,
                -0.007,
                Some(Element::H),
                epos[172],
                evel[172],
                None
            ),
            atom!(
                "H15Y",
                174,
                "POPC",
                3,
                1.008,
                -0.007,
                Some(Element::H),
                epos[173],
                evel[173],
                None
            ),
            atom!(
                "C316",
                175,
                "POPC",
                3,
                12.011,
                -0.081,
                Some(Element::C),
                epos[174],
                evel[174],
                None
            ),
            atom!(
                "H16X",
                176,
                "POPC",
                3,
                1.008,
                0.016,
                Some(Element::H),
                epos[175],
                evel[175],
                None
            ),
            atom!(
                "H16Y",
                177,
                "POPC",
                3,
                1.008,
                0.016,
                Some(Element::H),
                epos[176],
                evel[176],
                None
            ),
            atom!(
                "H16Z",
                178,
                "POPC",
                3,
                1.008,
                0.016,
                Some(Element::H),
                epos[177],
                evel[177],
                None
            ),
            // water
            atom!(
                "OW",
                179,
                "SOL",
                4,
                16.00,
                -0.834,
                Some(Element::O),
                epos[178],
                evel[178],
                None
            ),
            atom!(
                "HW1",
                180,
                "SOL",
                4,
                1.008,
                0.417,
                Some(Element::H),
                epos[179],
                evel[179],
                None
            ),
            atom!(
                "HW2",
                181,
                "SOL",
                4,
                1.008,
                0.417,
                Some(Element::H),
                epos[180],
                evel[180],
                None
            ),
            // cl ion
            atom!(
                "CL",
                182,
                "CL",
                5,
                35.45,
                -1.00,
                Some(Element::Cl),
                epos[181],
                evel[181],
                None
            ),
        ];

        for (a, e) in atoms.iter().zip(expected.iter()) {
            test_eq_atom(a, e);
        }

        let bonds = &topology.bonds;

        let mut expected_bonds = vec![
            bond!(0, 1),
            bond!(0, 2),
            bond!(0, 3),
            bond!(0, 4),
            bond!(4, 5),
            bond!(4, 6),
            bond!(4, 19),
            bond!(6, 7),
            bond!(6, 8),
            bond!(6, 9),
            bond!(9, 10),
            bond!(9, 11),
            bond!(9, 15),
            bond!(11, 12),
            bond!(11, 13),
            bond!(11, 14),
            bond!(15, 16),
            bond!(15, 17),
            bond!(15, 18),
            bond!(19, 20),
            bond!(19, 21),
            bond!(21, 22),
            bond!(21, 23),
            bond!(23, 24),
            bond!(23, 25),
            bond!(23, 41),
            bond!(25, 26),
            bond!(25, 27),
            bond!(25, 28),
            bond!(28, 29),
            bond!(28, 30),
            bond!(28, 31),
            bond!(31, 32),
            bond!(31, 33),
            bond!(31, 34),
            bond!(34, 35),
            bond!(34, 36),
            bond!(34, 37),
            bond!(37, 38),
            bond!(37, 39),
            bond!(37, 40),
            bond!(41, 42),
            bond!(41, 43),
            bond!(44, 45),
            bond!(44, 46),
            bond!(44, 47),
            bond!(44, 48),
            bond!(45, 49),
            bond!(45, 50),
            bond!(45, 60),
            bond!(46, 51),
            bond!(46, 52),
            bond!(46, 53),
            bond!(47, 54),
            bond!(47, 55),
            bond!(47, 56),
            bond!(48, 57),
            bond!(48, 58),
            bond!(48, 59),
            bond!(60, 61),
            bond!(60, 62),
            bond!(60, 66),
            bond!(63, 64),
            bond!(63, 65),
            bond!(63, 66),
            bond!(63, 67),
            bond!(67, 68),
            bond!(68, 69),
            bond!(68, 70),
            bond!(68, 71),
            bond!(71, 72),
            bond!(71, 73),
            bond!(71, 79),
            bond!(73, 74),
            bond!(74, 75),
            bond!(74, 76),
            bond!(76, 77),
            bond!(76, 78),
            bond!(76, 88),
            bond!(79, 80),
            bond!(79, 81),
            bond!(79, 82),
            bond!(82, 83),
            bond!(83, 84),
            bond!(83, 85),
            bond!(85, 86),
            bond!(85, 87),
            bond!(85, 135),
            bond!(88, 89),
            bond!(88, 90),
            bond!(88, 91),
            bond!(91, 92),
            bond!(91, 93),
            bond!(91, 94),
            bond!(94, 95),
            bond!(94, 96),
            bond!(94, 97),
            bond!(97, 98),
            bond!(97, 99),
            bond!(97, 100),
            bond!(100, 101),
            bond!(100, 102),
            bond!(100, 103),
            bond!(103, 104),
            bond!(103, 105),
            bond!(103, 106),
            bond!(106, 107),
            bond!(106, 108),
            bond!(108, 109),
            bond!(108, 110),
            bond!(110, 111),
            bond!(110, 112),
            bond!(110, 113),
            bond!(113, 114),
            bond!(113, 115),
            bond!(113, 116),
            bond!(116, 117),
            bond!(116, 118),
            bond!(116, 119),
            bond!(119, 120),
            bond!(119, 121),
            bond!(119, 122),
            bond!(122, 123),
            bond!(122, 124),
            bond!(122, 125),
            bond!(125, 126),
            bond!(125, 127),
            bond!(125, 128),
            bond!(128, 129),
            bond!(128, 130),
            bond!(128, 131),
            bond!(131, 132),
            bond!(131, 133),
            bond!(131, 134),
            bond!(135, 136),
            bond!(135, 137),
            bond!(135, 138),
            bond!(138, 139),
            bond!(138, 140),
            bond!(138, 141),
            bond!(141, 142),
            bond!(141, 143),
            bond!(141, 144),
            bond!(144, 145),
            bond!(144, 146),
            bond!(144, 147),
            bond!(147, 148),
            bond!(147, 149),
            bond!(147, 150),
            bond!(150, 151),
            bond!(150, 152),
            bond!(150, 153),
            bond!(153, 154),
            bond!(153, 155),
            bond!(153, 156),
            bond!(156, 157),
            bond!(156, 158),
            bond!(156, 159),
            bond!(159, 160),
            bond!(159, 161),
            bond!(159, 162),
            bond!(162, 163),
            bond!(162, 164),
            bond!(162, 165),
            bond!(165, 166),
            bond!(165, 167),
            bond!(165, 168),
            bond!(168, 169),
            bond!(168, 170),
            bond!(168, 171),
            bond!(171, 172),
            bond!(171, 173),
            bond!(171, 174),
            bond!(174, 175),
            bond!(174, 176),
            bond!(174, 177),
            bond!(178, 179),
            bond!(178, 180),
        ];

        if intermolecular {
            expected_bonds.extend(vec![bond!(15, 135), bond!(179, 181), bond!(149, 178)]);
        }

        assert_eq!(bonds.len(), expected_bonds.len());

        for e in expected_bonds.iter() {
            assert!(bonds.contains(e));
        }
    }

    #[test]
    fn small_aa_5() {
        let tpr = TprFile::parse("tests/test_files/small_aa_5.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 5.1.4"));
        assert_eq!(header.tpr_version, 103);
        assert_eq!(header.tpr_generation, 26);
        assert!(header.body_size.is_none());

        test_eq_small_aa(&tpr, false, GmxVersion::Gromacs5);
    }

    #[test]
    fn small_aa_2016() {
        let tpr = TprFile::parse("tests/test_files/small_aa_2016.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2016.4"));
        assert_eq!(header.tpr_version, 110);
        assert_eq!(header.tpr_generation, 26);
        assert!(header.body_size.is_none());

        test_eq_small_aa(&tpr, false, GmxVersion::Gromacs2016);
    }

    #[test]
    fn small_aa_2021() {
        let tpr = TprFile::parse("tests/test_files/small_aa_2021.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2021.4"));
        assert_eq!(header.tpr_version, 122);
        assert_eq!(header.tpr_generation, 28);
        assert_eq!(header.body_size.unwrap(), 70119);

        test_eq_small_aa(&tpr, false, GmxVersion::Gromacs2021);
    }

    #[test]
    fn small_aa_5_intermolecular() {
        let tpr = TprFile::parse("tests/test_files/small_aa_5_intermolecular.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 5.1.4"));
        assert_eq!(header.tpr_version, 103);
        assert_eq!(header.tpr_generation, 26);
        assert!(header.body_size.is_none());

        test_eq_small_aa(&tpr, true, GmxVersion::Gromacs5);
    }

    #[test]
    fn small_aa_2016_intermolecular() {
        let tpr = TprFile::parse("tests/test_files/small_aa_2016_intermolecular.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2016.4"));
        assert_eq!(header.tpr_version, 110);
        assert_eq!(header.tpr_generation, 26);
        assert!(header.body_size.is_none());

        test_eq_small_aa(&tpr, true, GmxVersion::Gromacs2016);
    }

    #[test]
    fn small_aa_2021_intermolecular() {
        let tpr = TprFile::parse("tests/test_files/small_aa_2021_intermolecular.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2021.4"));
        assert_eq!(header.tpr_version, 122);
        assert_eq!(header.tpr_generation, 28);
        assert_eq!(header.body_size.unwrap(), 70639);

        test_eq_small_aa(&tpr, true, GmxVersion::Gromacs2021);
    }

    #[test]
    fn large_5() {
        let tpr = TprFile::parse("tests/test_files/large_5.tpr").unwrap();

        assert_eq!(tpr.header.n_atoms, 7782);

        let first_atom = atom!(
            "BB",
            1,
            "LEU",
            1,
            72.0,
            1.0,
            None,
            Some([4.660, 1.341, 7.436]),
            Some([0.14950042963027954, -0.223041370511055, -0.2777465283870697]),
            None
        );
        let last_atom = atom!(
            "CL",
            7782,
            "ION",
            4926,
            72.0,
            -1.0,
            None,
            Some([2.753, 3.055, 1.049]),
            Some([
                -0.02956967242062092,
                0.05761498957872391,
                -0.43169528245925903
            ]),
            None
        );

        let first_bond = bond!(0, 1);
        let last_bond = bond!(3154, 3155);

        test_eq_atom(&tpr.topology.atoms[0], &first_atom);
        test_eq_atom(&tpr.topology.atoms[7781], &last_atom);

        assert!(&tpr.topology.bonds.contains(&first_bond));
        assert!(&tpr.topology.bonds.contains(&last_bond));
    }

    #[test]
    fn large_5_posres() {
        let tpr = TprFile::parse("tests/test_files/large_5_posres.tpr").unwrap();

        assert_eq!(tpr.header.n_atoms, 32443);

        let first_atom = atom!(
            "BB",
            1,
            "MET",
            1,
            72.0,
            1.0,
            None,
            Some([9.255, 5.545, 7.369]),
            Some([
                -0.09143112599849701,
                0.023182619363069534,
                0.42316368222236633
            ]),
            None
        );
        let last_atom = atom!(
            "CL",
            32443,
            "ION",
            21591,
            35.453,
            -1.0,
            None,
            Some([7.111, 9.099, 9.390]),
            Some([
                0.32411545515060425,
                -0.08399730175733566,
                0.1225026398897171
            ]),
            None
        );

        let first_bond = bond!(0, 2);
        let last_bond = bond!(12336, 12337);

        test_eq_atom(&tpr.topology.atoms[0], &first_atom);
        test_eq_atom(&tpr.topology.atoms[32442], &last_atom);

        assert!(&tpr.topology.bonds.contains(&first_bond));
        assert!(&tpr.topology.bonds.contains(&last_bond));
    }

    #[test]
    fn large_2021() {
        let tpr = TprFile::parse("tests/test_files/large_2021.tpr").unwrap();

        assert_eq!(tpr.header.n_atoms, 16949);

        let first_atom = atom!(
            "BB",
            1,
            "GLY",
            1,
            72.0,
            1.0,
            None,
            Some([5.788, 6.567, 7.798]),
            Some([
                0.18527933955192566,
                0.08655133843421936,
                -0.24502147734165192
            ]),
            None
        );
        let last_atom = atom!(
            "CL",
            16949,
            "ION",
            11209,
            35.453,
            -1.0,
            None,
            Some([4.271, 1.849, 10.175]),
            Some([
                -0.28008925914764404,
                0.6363704204559326,
                -0.030330466106534004
            ]),
            None
        );

        let first_bond = bond!(1, 2);
        let last_bond = bond!(6308, 6309);

        test_eq_atom(&tpr.topology.atoms[0], &first_atom);
        test_eq_atom(&tpr.topology.atoms[16948], &last_atom);

        assert!(&tpr.topology.bonds.contains(&first_bond));
        assert!(&tpr.topology.bonds.contains(&last_bond));
    }

    #[test]
    fn large_2021_aa() {
        let tpr = TprFile::parse("tests/test_files/large_2021_aa.tpr").unwrap();

        assert_eq!(tpr.header.n_atoms, 32817);

        let first_atom = atom!(
            "N",
            1,
            "SER",
            1,
            14.01,
            0.1849,
            Some(Element::N),
            Some([3.791, 4.267, 5.134]),
            Some([0.38751497864723206, -0.3888101875782013, 0.2692749500274658]),
            None
        );
        let last_atom = atom!(
            "CL",
            32817,
            "CL",
            5271,
            35.45,
            -1.0,
            Some(Element::Cl),
            Some([5.034, 2.0367, 6.9884]),
            Some([
                -0.03773900493979454,
                -0.1503848433494568,
                0.14367437362670898
            ]),
            None
        );

        let first_bond = bond!(0, 1);
        let last_bond = bond!(32785, 32787);

        test_eq_atom(&tpr.topology.atoms[0], &first_atom);
        test_eq_atom(&tpr.topology.atoms[32816], &last_atom);

        assert!(&tpr.topology.bonds.contains(&first_bond));
        assert!(&tpr.topology.bonds.contains(&last_bond));
    }

    #[test]
    fn large_2021_aa_posres() {
        let tpr = TprFile::parse("tests/test_files/large_2021_aa_posres.tpr").unwrap();

        assert_eq!(tpr.header.n_atoms, 34466);

        let first_atom = atom!(
            "N",
            1,
            "POPC",
            1,
            14.007,
            -0.6,
            Some(Element::N),
            Some([3.777, 3.390, 6.110]),
            Some([0.5271854400634766, -0.4534889757633209, 0.2048071026802063]),
            None
        );
        let last_atom = atom!(
            "H2",
            34466,
            "TIP3",
            5918,
            1.008,
            0.417,
            Some(Element::H),
            Some([5.742, 5.719, 2.157]),
            Some([
                -2.7277934551239014,
                0.46639126539230347,
                0.40979987382888794
            ]),
            None
        );

        let first_bond = bond!(0, 1);
        let last_bond = bond!(34463, 34465);

        test_eq_atom(&tpr.topology.atoms[0], &first_atom);
        test_eq_atom(&tpr.topology.atoms[34465], &last_atom);

        assert!(&tpr.topology.bonds.contains(&first_bond));
        assert!(&tpr.topology.bonds.contains(&last_bond));
    }

    #[test]
    fn double_2023() {
        let tpr = TprFile::parse("tests/test_files/double_2023.tpr").unwrap();

        assert_eq!(tpr.header.n_atoms, 16844);
        assert_eq!(tpr.header.gromacs_version, String::from("VERSION 2023.2"));
        assert_eq!(tpr.header.tpr_version, 129);
        assert_eq!(tpr.header.tpr_generation, 28);
        assert_eq!(tpr.header.body_size.unwrap(), 848223);

        let first_atom = atom!(
            "BB",
            1,
            "GLY",
            1,
            72.0,
            1.0,
            None,
            Some([9.497, 1.989, 7.498]),
            Some([
                -0.03727273095103525,
                0.07257158732040796,
                -0.00162505829335557
            ]),
            None
        );
        let last_atom = atom!(
            "CL",
            16844,
            "ION",
            11180,
            35.453,
            -1.0,
            None,
            Some([8.829, 11.186, 2.075]),
            Some([
                0.16692737861939136,
                0.16741210040248114,
                -0.27088445254642673
            ]),
            None
        );

        let first_bond = bond!(1, 2); // bond!(0, 1) is also present
        let last_bond = bond!(6203, 6204);

        test_eq_atom(&tpr.topology.atoms[0], &first_atom);
        test_eq_atom(&tpr.topology.atoms[16843], &last_atom);

        assert!(&tpr.topology.bonds.contains(&first_bond));
        assert!(&tpr.topology.bonds.contains(&last_bond));
    }

    fn test_eq_triclinic(tpr: &TprFile) {
        let simbox = tpr.simbox.as_ref().unwrap();

        assert_approx_eq!(f64, simbox.simbox[0][0], 5.29700, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[0][1], 0.00000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[0][2], 0.00000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox[1][0], 0.84445, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[1][1], 4.78912, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[1][2], 0.00000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox[2][0], 1.01785, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[2][1], -1.69043, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox[2][2], 2.22778, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_rel[0][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[0][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[0][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_rel[1][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[1][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[1][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_rel[2][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[2][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_rel[2][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_v[0][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[0][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[0][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_v[1][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[1][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[1][2], 0.0000, epsilon = 0.000001);

        assert_approx_eq!(f64, simbox.simbox_v[2][0], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[2][1], 0.0000, epsilon = 0.000001);
        assert_approx_eq!(f64, simbox.simbox_v[2][2], 0.0000, epsilon = 0.000001);

        assert_eq!(tpr.header.n_atoms, 50);

        let first_atom = atom!(
            "BB",
            1,
            "THR",
            1,
            72.0,
            1.0,
            None,
            Some([2.197, 0.567, 1.224]),
            Some([0.0, 0.0, 0.0]),
            None
        );
        let last_atom = atom!(
            "SC2",
            50,
            "LYS",
            21,
            54.0,
            1.0,
            None,
            Some([5.243, 2.953, 0.884]),
            Some([0.0, 0.0, 0.0]),
            None
        );

        let first_bond = bond!(0, 1);
        let last_bond = bond!(48, 49);

        test_eq_atom(&tpr.topology.atoms[0], &first_atom);
        test_eq_atom(&tpr.topology.atoms[49], &last_atom);

        assert!(&tpr.topology.bonds.contains(&first_bond));
        assert!(&tpr.topology.bonds.contains(&last_bond));
    }

    #[test]
    fn triclinic_2021() {
        let tpr = TprFile::parse("tests/test_files/triclinic_2021.tpr").unwrap();
        test_eq_triclinic(&tpr);
    }

    #[test]
    fn triclinic_5() {
        let tpr = TprFile::parse("tests/test_files/triclinic_5.tpr").unwrap();
        test_eq_triclinic(&tpr);
    }

    #[test]
    fn water_2021() {
        let tpr = TprFile::parse("tests/test_files/water_2021.tpr").unwrap();

        assert_eq!(tpr.topology.atoms.len(), 9);
        assert_eq!(tpr.topology.bonds.len(), 6);

        let expected_atom_names = ["OH2", "H1", "H2", "OH2", "H1", "H2", "OH2", "H1", "H2"];
        let expected_bonds = [(0, 1), (0, 2), (3, 4), (3, 5), (6, 7), (6, 8)];

        for (atom, expected) in tpr
            .topology
            .atoms
            .iter()
            .zip(expected_atom_names.into_iter())
        {
            assert_eq!(atom.atom_name, expected);
        }

        for (bond, expected) in tpr.topology.bonds.iter().zip(expected_bonds.into_iter()) {
            assert_eq!(bond.atom1, expected.0);
            assert_eq!(bond.atom2, expected.1);
        }
    }

    #[test]
    fn restrangles_2025() {
        let tpr = TprFile::parse("tests/test_files/restrangles_2025.tpr").unwrap();

        assert_eq!(tpr.topology.atoms.len(), 24404);
        assert_eq!(tpr.topology.bonds.len(), 9184);

        assert_eq!(
            tpr.topology.atoms.first().unwrap().atom_name,
            String::from("BB")
        );
        assert_eq!(
            tpr.topology.atoms.last().unwrap().atom_name,
            String::from("CL")
        );

        assert_eq!(tpr.topology.bonds.first().unwrap().atom1, 0);
        assert_eq!(tpr.topology.bonds.first().unwrap().atom2, 2);
        assert_eq!(tpr.topology.bonds.last().unwrap().atom1, 8495);
        assert_eq!(tpr.topology.bonds.last().unwrap().atom2, 8496);
    }
}

#[cfg(test)]
#[cfg(feature = "serde")]
mod tests_serde {
    use super::test_utilities::*;
    use minitpr::TprFile;
    use std::fs::read_to_string;

    #[test]
    fn to_yaml() {
        let tpr = TprFile::parse("tests/test_files/small_aa_2021.tpr").unwrap();

        let string = serde_yaml::to_string(&tpr).unwrap();
        let expected = read_to_string("tests/test_files/small_aa_2021.yaml").unwrap();

        assert_eq!(string, expected);
    }

    #[test]
    fn from_yaml() {
        let expected = TprFile::parse("tests/test_files/small_aa_2021.tpr").unwrap();
        let from_yaml: TprFile =
            serde_yaml::from_str(&read_to_string("tests/test_files/small_aa_2021.yaml").unwrap())
                .unwrap();

        for (a, e) in from_yaml
            .topology
            .atoms
            .iter()
            .zip(expected.topology.atoms.iter())
        {
            test_eq_atom(a, e);
        }

        assert_eq!(from_yaml.topology.bonds, expected.topology.bonds);
    }
}
