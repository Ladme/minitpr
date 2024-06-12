// Released under Apache License 2.0 / MIT License.
// Copyright (c) 2024 Ladislav Bartos

#[macro_use]
mod test_utilities {
    use float_cmp::assert_approx_eq;
    use minitpr::Atom;

    macro_rules! atom {
        ($atom_name:expr, $atom_number:expr, $residue_name:expr, $residue_number:expr, $mass:expr, $charge:expr, $element:expr) => {
            Atom {
                atom_name: $atom_name.to_owned(),
                atom_number: $atom_number,
                residue_name: $residue_name.to_owned(),
                residue_number: $residue_number,
                mass: $mass,
                charge: $charge,
                element: $element,
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

    pub(super) fn test_eq_atom(atom: &Atom, expected: &Atom) {
        assert_eq!(atom.atom_name, expected.atom_name);
        assert_eq!(atom.atom_number, expected.atom_number);
        assert_eq!(atom.residue_name, expected.residue_name);
        assert_eq!(atom.residue_number, expected.residue_number);
        assert_approx_eq!(f64, atom.mass, expected.mass, epsilon = 0.000001);
        assert_approx_eq!(f64, atom.charge, expected.charge, epsilon = 0.000001);
        assert_eq!(atom.element, expected.element);
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
        assert!(header.has_coordinates);
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

        let expected = vec![
            // Peptide
            atom!("BB", 1, "LEU", 1, 72.0, 1.0, None),
            atom!("SC1", 2, "LEU", 1, 54.0, 0.0, None),
            atom!("BB", 3, "SER", 2, 72.0, 0.0, None),
            atom!("SC1", 4, "SER", 2, 36.0, 0.0, None),
            atom!("BB", 5, "SER", 3, 72.0, 0.0, None),
            atom!("SC1", 6, "SER", 3, 36.0, 0.0, None),
            atom!("BB", 7, "LEU", 4, 72.0, 0.0, None),
            atom!("SC1", 8, "LEU", 4, 54.0, 0.0, None),
            atom!("BB", 9, "LEU", 5, 72.0, 0.0, None),
            atom!("SC1", 10, "LEU", 5, 54.0, 0.0, None),
            atom!("BB", 11, "SER", 6, 72.0, 0.0, None),
            atom!("SC1", 12, "SER", 6, 36.0, 0.0, None),
            atom!("BB", 13, "LEU", 7, 72.0, 0.0, None),
            atom!("SC1", 14, "LEU", 7, 54.0, 0.0, None),
            atom!("BB", 15, "LEU", 8, 72.0, 0.0, None),
            atom!("SC1", 16, "LEU", 8, 54.0, 0.0, None),
            atom!("BB", 17, "SER", 9, 72.0, 0.0, None),
            atom!("SC1", 18, "SER", 9, 36.0, 0.0, None),
            atom!("BB", 19, "SER", 10, 72.0, 0.0, None),
            atom!("SC1", 20, "SER", 10, 36.0, 0.0, None),
            atom!("BB", 21, "LEU", 11, 72.0, 0.0, None),
            atom!("SC1", 22, "LEU", 11, 54.0, 0.0, None),
            atom!("BB", 23, "LEU", 12, 72.0, 0.0, None),
            atom!("SC1", 24, "LEU", 12, 54.0, 0.0, None),
            atom!("BB", 25, "SER", 13, 72.0, 0.0, None),
            atom!("SC1", 26, "SER", 13, 36.0, 0.0, None),
            atom!("BB", 27, "LEU", 14, 72.0, 0.0, None),
            atom!("SC1", 28, "LEU", 14, 54.0, 0.0, None),
            atom!("BB", 29, "LEU", 15, 72.0, 0.0, None),
            atom!("SC1", 30, "LEU", 15, 54.0, 0.0, None),
            atom!("BB", 31, "SER", 16, 72.0, 0.0, None),
            atom!("SC1", 32, "SER", 16, 36.0, 0.0, None),
            atom!("BB", 33, "SER", 17, 72.0, 0.0, None),
            atom!("SC1", 34, "SER", 17, 36.0, 0.0, None),
            atom!("BB", 35, "LEU", 18, 72.0, 0.0, None),
            atom!("SC1", 36, "LEU", 18, 54.0, 0.0, None),
            atom!("BB", 37, "LEU", 19, 72.0, 0.0, None),
            atom!("SC1", 38, "LEU", 19, 54.0, 0.0, None),
            atom!("BB", 39, "SER", 20, 72.0, 0.0, None),
            atom!("SC1", 40, "SER", 20, 36.0, 0.0, None),
            atom!("BB", 41, "LEU", 21, 72.0, 0.0, None),
            atom!("SC1", 42, "LEU", 21, 54.0, 0.0, None),
            // POPC #1
            atom!("NC3", 43, "POPC", 22, 72.0, 1.0, None),
            atom!("PO4", 44, "POPC", 22, 72.0, -1.0, None),
            atom!("GL1", 45, "POPC", 22, 54.0, 0.0, None),
            atom!("GL2", 46, "POPC", 22, 72.0, 0.0, None),
            atom!("C1A", 47, "POPC", 22, 72.0, 0.0, None),
            atom!("D2A", 48, "POPC", 22, 72.0, 0.0, None),
            atom!("C3A", 49, "POPC", 22, 72.0, 0.0, None),
            atom!("C4A", 50, "POPC", 22, 72.0, 0.0, None),
            atom!("C1B", 51, "POPC", 22, 72.0, 0.0, None),
            atom!("C2B", 52, "POPC", 22, 72.0, 0.0, None),
            atom!("C3B", 53, "POPC", 22, 72.0, 0.0, None),
            atom!("C4B", 54, "POPC", 22, 72.0, 0.0, None),
            // POPC #2
            atom!("NC3", 55, "POPC", 23, 72.0, 1.0, None),
            atom!("PO4", 56, "POPC", 23, 72.0, -1.0, None),
            atom!("GL1", 57, "POPC", 23, 54.0, 0.0, None),
            atom!("GL2", 58, "POPC", 23, 72.0, 0.0, None),
            atom!("C1A", 59, "POPC", 23, 72.0, 0.0, None),
            atom!("D2A", 60, "POPC", 23, 72.0, 0.0, None),
            atom!("C3A", 61, "POPC", 23, 72.0, 0.0, None),
            atom!("C4A", 62, "POPC", 23, 72.0, 0.0, None),
            atom!("C1B", 63, "POPC", 23, 72.0, 0.0, None),
            atom!("C2B", 64, "POPC", 23, 72.0, 0.0, None),
            atom!("C3B", 65, "POPC", 23, 72.0, 0.0, None),
            atom!("C4B", 66, "POPC", 23, 72.0, 0.0, None),
            // water (10x)
            atom!("W", 67, "W", 24, 72.0, 0.0, None),
            atom!("W", 68, "W", 25, 72.0, 0.0, None),
            atom!("W", 69, "W", 26, 72.0, 0.0, None),
            atom!("W", 70, "W", 27, 72.0, 0.0, None),
            atom!("W", 71, "W", 28, 72.0, 0.0, None),
            atom!("W", 72, "W", 29, 72.0, 0.0, None),
            atom!("W", 73, "W", 30, 72.0, 0.0, None),
            atom!("W", 74, "W", 31, 72.0, 0.0, None),
            atom!("W", 75, "W", 32, 72.0, 0.0, None),
            atom!("W", 76, "W", 33, 72.0, 0.0, None),
            // CL ion
            atom!("CL-", 77, "ION", 34, 35.453, -1.0, None),
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
        match TprFile::parse("tests/test_files/empty.tpr") {
            Ok(_) => panic!("Parsing should have failed."),
            Err(_) => (),
        }
    }

    fn test_eq_small_aa(tpr: &TprFile, intermolecular: bool) {
        let header = &tpr.header;

        assert_eq!(header.precision, Precision::Single);
        assert_eq!(header.file_tag, "release");
        assert_eq!(header.n_atoms, 182);
        assert_eq!(header.n_coupling_groups, 3);
        assert_eq!(header.fep_state, 0);
        assert_approx_eq!(f64, header.lambda, 0.0);
        assert!(header.has_input_record);
        assert!(header.has_topology);
        assert!(header.has_coordinates);
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

        let expected = vec![
            // dipeptide
            atom!("N", 1, "LEU", 1, 14.01, 0.101, Some(Element::N)),
            atom!("H1", 2, "LEU", 1, 1.008, 0.2148, Some(Element::H)),
            atom!("H2", 3, "LEU", 1, 1.008, 0.2148, Some(Element::H)),
            atom!("H3", 4, "LEU", 1, 1.008, 0.2148, Some(Element::H)),
            atom!("CA", 5, "LEU", 1, 12.01, 0.0104, Some(Element::C)),
            atom!("HA", 6, "LEU", 1, 1.008, 0.1053, Some(Element::H)),
            atom!("CB", 7, "LEU", 1, 12.01, -0.0244, Some(Element::C)),
            atom!("HB1", 8, "LEU", 1, 1.008, 0.0256, Some(Element::H)),
            atom!("HB2", 9, "LEU", 1, 1.008, 0.0256, Some(Element::H)),
            atom!("CG", 10, "LEU", 1, 12.01, 0.3421, Some(Element::C)),
            atom!("HG", 11, "LEU", 1, 1.008, -0.038, Some(Element::H)),
            atom!("CD1", 12, "LEU", 1, 12.01, -0.4106, Some(Element::C)),
            atom!("HD11", 13, "LEU", 1, 1.008, 0.098, Some(Element::H)),
            atom!("HD12", 14, "LEU", 1, 1.008, 0.098, Some(Element::H)),
            atom!("HD13", 15, "LEU", 1, 1.008, 0.098, Some(Element::H)),
            atom!("CD2", 16, "LEU", 1, 12.01, -0.4104, Some(Element::C)),
            atom!("HD21", 17, "LEU", 1, 1.008, 0.098, Some(Element::H)),
            atom!("HD22", 18, "LEU", 1, 1.008, 0.098, Some(Element::H)),
            atom!("HD23", 19, "LEU", 1, 1.008, 0.098, Some(Element::H)),
            atom!("C", 20, "LEU", 1, 12.01, 0.6123, Some(Element::C)),
            atom!("O", 21, "LEU", 1, 16.00, -0.5713, Some(Element::O)),
            atom!("N", 22, "LYS", 2, 14.01, -0.3481, Some(Element::N)),
            atom!("H", 23, "LYS", 2, 1.008, 0.2764, Some(Element::H)),
            atom!("CA", 24, "LYS", 2, 12.01, -0.2903, Some(Element::C)),
            atom!("HA", 25, "LYS", 2, 1.008, 0.1438, Some(Element::H)),
            atom!("CB", 26, "LYS", 2, 12.01, -0.0538, Some(Element::C)),
            atom!("HB1", 27, "LYS", 2, 1.008, 0.0482, Some(Element::H)),
            atom!("HB2", 28, "LYS", 2, 1.008, 0.0482, Some(Element::H)),
            atom!("CG", 29, "LYS", 2, 12.01, 0.0227, Some(Element::C)),
            atom!("HG1", 30, "LYS", 2, 1.008, 0.0134, Some(Element::H)),
            atom!("HG2", 31, "LYS", 2, 1.008, 0.0134, Some(Element::H)),
            atom!("CD", 32, "LYS", 2, 12.01, -0.0392, Some(Element::C)),
            atom!("HD1", 33, "LYS", 2, 1.008, 0.0611, Some(Element::H)),
            atom!("HD2", 34, "LYS", 2, 1.008, 0.0611, Some(Element::H)),
            atom!("CE", 35, "LYS", 2, 12.01, -0.0176, Some(Element::C)),
            atom!("HE1", 36, "LYS", 2, 1.008, 0.1121, Some(Element::H)),
            atom!("HE2", 37, "LYS", 2, 1.008, 0.1121, Some(Element::H)),
            atom!("NZ", 38, "LYS", 2, 14.01, -0.3741, Some(Element::N)),
            atom!("HZ1", 39, "LYS", 2, 1.008, 0.3374, Some(Element::H)),
            atom!("HZ2", 40, "LYS", 2, 1.008, 0.3374, Some(Element::H)),
            atom!("HZ3", 41, "LYS", 2, 1.008, 0.3374, Some(Element::H)),
            atom!("C", 42, "LYS", 2, 12.01, 0.8488, Some(Element::C)),
            atom!("OC1", 43, "LYS", 2, 16.00, -0.8252, Some(Element::O)),
            atom!("OC2", 44, "LYS", 2, 16.00, -0.8252, Some(Element::O)),
            // POPC
            atom!("N", 45, "POPC", 3, 14.007, 0.20, Some(Element::N)),
            atom!("C12", 46, "POPC", 3, 12.011, -0.20, Some(Element::C)),
            atom!("C13", 47, "POPC", 3, 12.011, -0.38, Some(Element::C)),
            atom!("C14", 48, "POPC", 3, 12.011, -0.38, Some(Element::C)),
            atom!("C15", 49, "POPC", 3, 12.011, -0.38, Some(Element::C)),
            atom!("H12A", 50, "POPC", 3, 1.008, 0.09, Some(Element::H)),
            atom!("H12B", 51, "POPC", 3, 1.008, 0.09, Some(Element::H)),
            atom!("H13A", 52, "POPC", 3, 1.008, 0.19, Some(Element::H)),
            atom!("H13B", 53, "POPC", 3, 1.008, 0.19, Some(Element::H)),
            atom!("H13C", 54, "POPC", 3, 1.008, 0.19, Some(Element::H)),
            atom!("H14A", 55, "POPC", 3, 1.008, 0.19, Some(Element::H)),
            atom!("H14B", 56, "POPC", 3, 1.008, 0.19, Some(Element::H)),
            atom!("H14C", 57, "POPC", 3, 1.008, 0.19, Some(Element::H)),
            atom!("H15A", 58, "POPC", 3, 1.008, 0.19, Some(Element::H)),
            atom!("H15B", 59, "POPC", 3, 1.008, 0.19, Some(Element::H)),
            atom!("H15C", 60, "POPC", 3, 1.008, 0.19, Some(Element::H)),
            atom!("C11", 61, "POPC", 3, 12.011, 0.17, Some(Element::C)),
            atom!("H11A", 62, "POPC", 3, 1.008, 0.03, Some(Element::H)),
            atom!("H11B", 63, "POPC", 3, 1.008, 0.03, Some(Element::H)),
            atom!("P", 64, "POPC", 3, 30.974, 1.58, Some(Element::P)),
            atom!("O13", 65, "POPC", 3, 15.9994, -0.86, Some(Element::O)),
            atom!("O14", 66, "POPC", 3, 15.9994, -0.86, Some(Element::O)),
            atom!("O12", 67, "POPC", 3, 15.9994, -0.49, Some(Element::O)),
            atom!("O11", 68, "POPC", 3, 15.9994, -0.49, Some(Element::O)),
            atom!("C1", 69, "POPC", 3, 12.011, -0.11, Some(Element::C)),
            atom!("HA", 70, "POPC", 3, 1.008, 0.07, Some(Element::H)),
            atom!("HB", 71, "POPC", 3, 1.008, 0.07, Some(Element::H)),
            atom!("C2", 72, "POPC", 3, 12.011, 0.48, Some(Element::C)),
            atom!("HS", 73, "POPC", 3, 1.008, 0.04, Some(Element::H)),
            atom!("O21", 74, "POPC", 3, 15.9994, -0.47, Some(Element::O)),
            atom!("C21", 75, "POPC", 3, 12.011, 0.79, Some(Element::C)),
            atom!("O22", 76, "POPC", 3, 15.9994, -0.65, Some(Element::O)),
            atom!("C22", 77, "POPC", 3, 12.011, -0.06, Some(Element::C)),
            atom!("H2R", 78, "POPC", 3, 1.008, 0.03, Some(Element::H)),
            atom!("H2S", 79, "POPC", 3, 1.008, 0.03, Some(Element::H)),
            atom!("C3", 80, "POPC", 3, 12.011, 0.13, Some(Element::C)),
            atom!("HX", 81, "POPC", 3, 1.008, 0.06, Some(Element::H)),
            atom!("HY", 82, "POPC", 3, 1.008, 0.06, Some(Element::H)),
            atom!("O31", 83, "POPC", 3, 15.9994, -0.47, Some(Element::O)),
            atom!("C31", 84, "POPC", 3, 12.011, 0.79, Some(Element::C)),
            atom!("O32", 85, "POPC", 3, 15.9994, -0.65, Some(Element::O)),
            atom!("C32", 86, "POPC", 3, 12.011, -0.06, Some(Element::C)),
            atom!("H2X", 87, "POPC", 3, 1.008, 0.03, Some(Element::H)),
            atom!("H2Y", 88, "POPC", 3, 1.008, 0.03, Some(Element::H)),
            atom!("C23", 89, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H3R", 90, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H3S", 91, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C24", 92, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H4R", 93, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H4S", 94, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C25", 95, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H5R", 96, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H5S", 97, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C26", 98, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H6R", 99, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H6S", 100, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C27", 101, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H7R", 102, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H7S", 103, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C28", 104, "POPC", 3, 12.011, 0.03, Some(Element::C)),
            atom!("H8R", 105, "POPC", 3, 1.008, 0.03, Some(Element::H)),
            atom!("H8S", 106, "POPC", 3, 1.008, 0.03, Some(Element::H)),
            atom!("C29", 107, "POPC", 3, 12.011, -0.20, Some(Element::C)),
            atom!("H91", 108, "POPC", 3, 1.008, 0.11, Some(Element::H)),
            atom!("C210", 109, "POPC", 3, 12.011, -0.20, Some(Element::C)),
            atom!("H101", 110, "POPC", 3, 1.008, 0.11, Some(Element::H)),
            atom!("C211", 111, "POPC", 3, 12.011, 0.03, Some(Element::C)),
            atom!("H11R", 112, "POPC", 3, 1.008, 0.03, Some(Element::H)),
            atom!("H11S", 113, "POPC", 3, 1.008, 0.03, Some(Element::H)),
            atom!("C212", 114, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H12R", 115, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H12S", 116, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C213", 117, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H13R", 118, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H13S", 119, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C214", 120, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H14R", 121, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H14S", 122, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C215", 123, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H15R", 124, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H15S", 125, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C216", 126, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H16R", 127, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H16S", 128, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C217", 129, "POPC", 3, 12.011, 0.047, Some(Element::C)),
            atom!("H17R", 130, "POPC", 3, 1.008, -0.007, Some(Element::H)),
            atom!("H17S", 131, "POPC", 3, 1.008, -0.007, Some(Element::H)),
            atom!("C218", 132, "POPC", 3, 12.011, -0.081, Some(Element::C)),
            atom!("H18R", 133, "POPC", 3, 1.008, 0.016, Some(Element::H)),
            atom!("H18S", 134, "POPC", 3, 1.008, 0.016, Some(Element::H)),
            atom!("H18T", 135, "POPC", 3, 1.008, 0.016, Some(Element::H)),
            atom!("C33", 136, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H3X", 137, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H3Y", 138, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C34", 139, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H4X", 140, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H4Y", 141, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C35", 142, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H5X", 143, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H5Y", 144, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C36", 145, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H6X", 146, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H6Y", 147, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C37", 148, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H7X", 149, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H7Y", 150, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C38", 151, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H8X", 152, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H8Y", 153, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C39", 154, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H9X", 155, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H9Y", 156, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C310", 157, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H10X", 158, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H10Y", 159, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C311", 160, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H11X", 161, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H11Y", 162, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C312", 163, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H12X", 164, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H12Y", 165, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C313", 166, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H13X", 167, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H13Y", 168, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C314", 169, "POPC", 3, 12.011, 0.00, Some(Element::C)),
            atom!("H14X", 170, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("H14Y", 171, "POPC", 3, 1.008, 0.00, Some(Element::H)),
            atom!("C315", 172, "POPC", 3, 12.011, 0.047, Some(Element::C)),
            atom!("H15X", 173, "POPC", 3, 1.008, -0.007, Some(Element::H)),
            atom!("H15Y", 174, "POPC", 3, 1.008, -0.007, Some(Element::H)),
            atom!("C316", 175, "POPC", 3, 12.011, -0.081, Some(Element::C)),
            atom!("H16X", 176, "POPC", 3, 1.008, 0.016, Some(Element::H)),
            atom!("H16Y", 177, "POPC", 3, 1.008, 0.016, Some(Element::H)),
            atom!("H16Z", 178, "POPC", 3, 1.008, 0.016, Some(Element::H)),
            // water
            atom!("OW", 179, "SOL", 4, 16.00, -0.834, Some(Element::O)),
            atom!("HW1", 180, "SOL", 4, 1.008, 0.417, Some(Element::H)),
            atom!("HW2", 181, "SOL", 4, 1.008, 0.417, Some(Element::H)),
            // cl ion
            atom!("CL", 182, "CL", 5, 35.45, -1.00, Some(Element::Cl)),
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

        test_eq_small_aa(&tpr, false);
    }

    #[test]
    fn small_aa_2016() {
        let tpr = TprFile::parse("tests/test_files/small_aa_2016.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2016.4"));
        assert_eq!(header.tpr_version, 110);
        assert_eq!(header.tpr_generation, 26);
        assert!(header.body_size.is_none());

        test_eq_small_aa(&tpr, false);
    }

    #[test]
    fn small_aa_2021() {
        let tpr = TprFile::parse("tests/test_files/small_aa_2021.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2021.4"));
        assert_eq!(header.tpr_version, 122);
        assert_eq!(header.tpr_generation, 28);
        assert_eq!(header.body_size.unwrap(), 70119);

        test_eq_small_aa(&tpr, false);
    }

    #[test]
    fn small_aa_5_intermolecular() {
        let tpr = TprFile::parse("tests/test_files/small_aa_5_intermolecular.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 5.1.4"));
        assert_eq!(header.tpr_version, 103);
        assert_eq!(header.tpr_generation, 26);
        assert!(header.body_size.is_none());

        test_eq_small_aa(&tpr, true);
    }

    #[test]
    fn small_aa_2016_intermolecular() {
        let tpr = TprFile::parse("tests/test_files/small_aa_2016_intermolecular.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2016.4"));
        assert_eq!(header.tpr_version, 110);
        assert_eq!(header.tpr_generation, 26);
        assert!(header.body_size.is_none());

        test_eq_small_aa(&tpr, true);
    }

    #[test]
    fn small_aa_2021_intermolecular() {
        let tpr = TprFile::parse("tests/test_files/small_aa_2021_intermolecular.tpr").unwrap();

        let header = &tpr.header;
        assert_eq!(header.gromacs_version, String::from("VERSION 2021.4"));
        assert_eq!(header.tpr_version, 122);
        assert_eq!(header.tpr_generation, 28);
        assert_eq!(header.body_size.unwrap(), 70639);

        test_eq_small_aa(&tpr, true);
    }

    #[test]
    fn large_5() {
        let tpr = TprFile::parse("tests/test_files/large_5.tpr").unwrap();

        assert_eq!(tpr.header.n_atoms, 7782);

        let first_atom = atom!("BB", 1, "LEU", 1, 72.0, 1.0, None);
        let last_atom = atom!("CL", 7782, "ION", 4926, 72.0, -1.0, None);

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

        let first_atom = atom!("BB", 1, "MET", 1, 72.0, 1.0, None);
        let last_atom = atom!("CL", 32443, "ION", 21591, 35.453, -1.0, None);

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

        let first_atom = atom!("BB", 1, "GLY", 1, 72.0, 1.0, None);
        let last_atom = atom!("CL", 16949, "ION", 11209, 35.453, -1.0, None);

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

        let first_atom = atom!("N", 1, "SER", 1, 14.01, 0.1849, Some(Element::N));
        let last_atom = atom!("CL", 32817, "CL", 5271, 35.45, -1.0, Some(Element::Cl));

        let first_bond = bond!(0, 1);
        let last_bond = bond!(17511, 17514);

        test_eq_atom(&tpr.topology.atoms[0], &first_atom);
        test_eq_atom(&tpr.topology.atoms[32816], &last_atom);

        assert!(&tpr.topology.bonds.contains(&first_bond));
        assert!(&tpr.topology.bonds.contains(&last_bond));
    }

    #[test]
    fn large_2021_aa_posres() {
        let tpr = TprFile::parse("tests/test_files/large_2021_aa_posres.tpr").unwrap();

        assert_eq!(tpr.header.n_atoms, 34466);

        let first_atom = atom!("N", 1, "POPC", 1, 14.007, -0.6, Some(Element::N));
        let last_atom = atom!("H2", 34466, "TIP3", 5918, 1.008, 0.417, Some(Element::H));

        let first_bond = bond!(0, 1);
        let last_bond = bond!(17148, 17151);

        test_eq_atom(&tpr.topology.atoms[0], &first_atom);
        test_eq_atom(&tpr.topology.atoms[34465], &last_atom);

        assert!(&tpr.topology.bonds.contains(&first_bond));
        assert!(&tpr.topology.bonds.contains(&last_bond));
    }

    #[test]
    fn double_2023() {
        let tpr = TprFile::parse("tests/test_files/double_2023.tpr").unwrap();

        assert_eq!(tpr.header.n_atoms, 16844);

        let first_atom = atom!("BB", 1, "GLY", 1, 72.0, 1.0, None);
        let last_atom = atom!("CL", 16844, "ION", 11180, 35.453, -1.0, None);

        let first_bond = bond!(1, 2);
        let last_bond = bond!(6203, 6204);

        test_eq_atom(&tpr.topology.atoms[0], &first_atom);
        test_eq_atom(&tpr.topology.atoms[16843], &last_atom);

        assert!(&tpr.topology.bonds.contains(&first_bond));
        assert!(&tpr.topology.bonds.contains(&last_bond));
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
        let tpr = TprFile::parse("tests/test_files/small_cg_2021.tpr").unwrap();

        let string = serde_yaml::to_string(&tpr).unwrap();
        let expected = read_to_string("tests/test_files/small_cg_2021.yaml").unwrap();

        assert_eq!(string, expected);
    }

    #[test]
    fn from_yaml() {
        let expected = TprFile::parse("tests/test_files/small_cg_2021.tpr").unwrap();
        let from_yaml: TprFile =
            serde_yaml::from_str(&read_to_string("tests/test_files/small_cg_2021.yaml").unwrap())
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
