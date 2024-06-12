# minitpr: Library for Reading Gromacs TPR Files


`minitpr` is a pure Rust library designed for parsing Gromacs tpr files, focusing on extracting system topology and structure.

## Capabilities and Limitations
- Supports parsing of tpr files from version 103 onwards (Gromacs 5.1 and later).
- Extracts system topology and structure: atoms, their basic properties (including positions, velocities, and forces), and bonds between atoms (including intermolecular bonds).
- Does **not** support parsing of force-field and simulation parameters, nor does it offer capabilities to write tpr files.

## Usage
To include `minitpr` in your project, add it as a dependency using Cargo:
```shell
cargo add minitpr
```

### Parsing a tpr file
Example usage to parse a tpr file and handle potential errors:
```rust
use minitpr::TprFile;

fn main() {
    let tpr = match TprFile::parse("topol.tpr") {
        Ok(file) => file,
        Err(error) => {
            eprintln!("{}", error);
            return;
        },
    };

    // now you can work with the `tpr` object, accessing its properties and data
    // for instance, to iterate through the atoms of the molecular system, use:
    for atom in tpr.topology.atoms.iter() {
        // perform some operation with the atom
    }
}
```

### Data Structures
The `TprFile` structure encapsulates the following information:

- Header: Metadata about the tpr file (see `TprHeader` structure).
- Molecular System Name: The name of the simulated system.
- Simulation Box Dimensions: Available within the `SimBox` structure if present.
- System Topology: Topology of the molecular system containing atoms and bonds (see `TprTopology` structure).

Each atom (see `Atom`) represented in the system topology includes:
- Atom name.
- Sequential atom number, starting from 1.
- Residue name.
- Sequential residue number, starting from 1.
- Mass.
- Charge.
- Element (`None` if unidentifiable).
- Position (`None` if not present).
- Velocity (`None` if not present).
- Force (`None` if not present).

## Features
### Serialization/Deserialization
Enable (de)serialization support for `TprFile` with `serde` by adding the feature flag during installation:
```shell
cargo add minitpr --features serde
```

## License
`minitpr` is open-sourced under either the [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0) or the [MIT License](https://opensource.org/license/MIT) at your option.

## Disclaimer
Due to the huge number of various force field and simulation parameters and Gromacs options, it is very difficult to *comprehensively* test the `minitpr` library.Your contributions in the form of tests, tpr files with unusual parameters, or new functionality are very welcome.

If the library is unable to parse your tpr file, but you believe it should be able to, please open a [GitHub issue](https://github.com/Ladme/minitpr/issues) and **upload your tpr file**.
