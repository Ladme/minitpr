# Changelog for the `minitpr` crate

## Version 0.2.2
- **BUG FIX**: Fixed bug where bonds were not loaded for some water models (TIP3P and similar). The SETTLE interaction is now properly translated into bonds.

## Version 0.2.0
- **BREAKING CHANGE:** `TprHeader` field `has_coordinates` has been renamed to `has_positions`.
- `minitpr` can now parse intermolecular bonds. Intermolecular bonds are added to the end of the `bonds` vector in `TprFile::TprTopology`.
- `minitpr` can now read positions, velocities, and forces of atoms.