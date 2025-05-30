## Overview

The `esl_geometry_m` module provides the `geometry_t` derived type, which serves as a comprehensive container for all geometric information defining a physical system for an electronic structure calculation. This includes details about the different atomic species present, the total number of atoms, their Cartesian coordinates (in Bohr units), and the unit cell vectors that define the simulation cell's size, shape, and periodicity. The module offers functionalities to initialize this geometric data (primarily by parsing an FDF input file), to calculate derived quantities like cell volume and total valence electronic charge, to provide a formatted summary of the geometry, and to ensure proper deallocation of its components.

## Key Components

- **Module:** `esl_geometry_m`
    - **Description:** Defines and manages the `geometry_t` derived type, which encapsulates all geometric aspects of the simulated system.
- **Type:** `geometry_t`
    - **Description:** A derived data type that holds structural information about the system, including atomic species, atomic coordinates, and unit cell parameters.
    - **Components:**
        - `n_species (integer(ip))`: The number of unique atomic species in the system. Default: `0`.
        - `species (type(species_t), allocatable :: species(:))`: An allocatable array of `species_t` objects (where `species_t` is defined in `esl_species_m`). Each element describes a unique atomic species (e.g., its label, pseudopotential, valence charge).
        - `n_atoms (integer(ip))`: The total number of atoms in the simulation cell. Default: `0`.
        - `xyz (real(dp), allocatable :: xyz(:,:))`: A 2D allocatable array storing the Cartesian coordinates of the atoms. Dimensions are `(3, n_atoms)`, with coordinates in Bohr.
        - `species_idx (integer(ip), allocatable :: species_idx(:))`: A 1D allocatable array that maps each atom (index from 1 to `n_atoms`) to its corresponding species type (index in the `species` array).
        - `cell (real(dp) :: cell(3,3))`: A 3x3 matrix where the columns represent the lattice vectors **a**, **b**, **c**. If the system is molecular (non-periodic), this might be set to define a bounding box. Default: all zeros.
        - `icell (real(dp) :: icell(3,3))`: The inverse of the `cell` matrix, useful for converting Cartesian coordinates to fractional coordinates. Default: all zeros.
        - `vol (real(dp))`: Stores the volume of the unit cell in Bohr^3. Default: `0.0_dp`.
    - **Procedures (public, type-bound):**
        - `init`: Initializes the geometry by reading data from an FDF input file.
        - `summary`: Prints a YAML-formatted summary of the geometry.
        - `volume`: Calculates the cell volume (and stores it in `this%vol`).
        - `electronic_charge`: Calculates the total number of valence electrons.
    - **Finalizer (type-bound):**
        - `finalizer`: Handles deallocation of allocatable arrays (`species`, `xyz`, `species_idx`).

### Type-Bound Procedures

#### `init(this)`
- **Description:** Initializes the `geometry_t` object by parsing data from an FDF (Flexible Data Format) input file. This comprehensive routine:
    1.  Reads unit cell parameters. It checks for a `Cubic` FDF key to define a cubic cell directly.
    2.  Reads atomic species information from the `Species` FDF block. Each entry typically defines a species label and associated properties (like pseudopotential file, mass, valence charge), which are handled by `this%species(is)%init()`. It errors if no species are defined.
    3.  Reads atomic coordinates from the `Coordinates` FDF block. It determines `n_atoms` (either by counting lines or from the `NumberOfAtoms` FDF key) and allocates `xyz` and `species_idx`. For each atom, it reads its x, y, z coordinates and a species label, then maps the label to an index in the `this%species` array. It errors if coordinates are missing or if an undefined species label is used.
    4.  If the cell was not explicitly defined (e.g., via `Cubic`), it constructs a default orthorhombic cell that tightly encloses all atoms, plus a vacuum padding of approximately 30 Bohr in each Cartesian direction. This is suitable for isolated molecule calculations.
    5.  Calculates the cell volume by calling `this%volume()` and stores it. It raises an error if the calculated volume is negative.
    6.  Computes the inverse cell matrix `this%icell` using `matr3inv(this%cell)`.
- **Arguments:** `this (class(geometry_t), intent(inout))`: The `geometry_t` object to be initialized.

#### `summary(this)`
- **Description:** Prints a human-readable summary of the geometry to standard output using YAML formatting. The summary includes:
    -   The cell vectors (`this%cell`).
    -   A list of atomic coordinates, showing the element label and its X, Y, Z coordinates for each atom.
    -   Information for each unique species (by calling `this%species(ia)%summary()`).
    -   The total number of valence electrons (obtained from `this%electronic_charge()`).
    -   The cell volume (obtained from `this%volume()`).
- **Arguments:** `this (class(geometry_t), intent(in))`: The `geometry_t` object to summarize.

#### `volume(this) result(vol_val)`
- **Description:** Calculates the volume of the unit cell defined by the `this%cell` lattice vectors using the scalar triple product: `vol = a . (b x c)`. The absolute value of this result is stored in `this%vol` and also returned by the function.
- **Arguments:** `this (class(geometry_t), intent(inout))`: The `geometry_t` object. `this%vol` is updated.
- **Returns:** `vol_val (real(dp))`: The absolute volume of the unit cell in Bohr^3.

#### `electronic_charge(this) result(charge)`
- **Description:** Calculates the total number of valence electrons in the system. It iterates through all atoms, identifies their species using `this%species_idx`, and sums the valence charge (`q` member of `species_t`) for each species.
- **Arguments:** `this (class(geometry_t), intent(in))`: The `geometry_t` object. *(Note: Source code has `intent(inout)`, but `in` should suffice as `this` is not modified)*.
- **Returns:** `charge (real(dp))`: The total valence electronic charge.

#### `finalizer(this)`
- **Description:** This is a type-bound `final` subroutine automatically invoked when a `geometry_t` object is deallocated. It ensures that the allocatable arrays `this%species`, `this%xyz`, and `this%species_idx` are deallocated if they were previously allocated, preventing memory leaks.
- **Arguments:** `this (type(geometry_t))`: The `geometry_t` object being finalized.

## Important Variables/Constants
- All components of `geometry_t` (e.g., `n_atoms`, `xyz`, `species`, `cell`) are critical for defining the system's physical structure.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage:
!
! use esl_geometry_m
! type(geometry_t) :: molecule_geom
! real(dp) :: total_valence_charge, cell_vol
!
! ! Initialize geometry from FDF input file (not shown: FDF setup)
! call molecule_geom%init()
!
! ! Get some information
! total_valence_charge = molecule_geom%electronic_charge()
! cell_vol = molecule_geom%volume() ! Also updates molecule_geom%vol
!
! ! Print summary
! call molecule_geom%summary()
!
! ! molecule_geom will be finalized automatically when it goes out of scope.
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec`: Provides precision kinds `dp` (double precision) and `ip` (integer).
    - `fdf`: This module is heavily relied upon by the `init` procedure to read all geometry-related information (cell parameters, species definitions, atomic coordinates) from an FDF-formatted input file.
    - `esl_numeric_m, only : matr3inv`: The `matr3inv` function is used to compute the inverse of the cell matrix during initialization.
    - `esl_species_m`: Provides the `species_t` derived type, which is used as an array component within `geometry_t` to store detailed information about each unique atomic species.
    - `esl_message_m`: Used for error reporting (e.g., `message_error` if input is missing or inconsistent).
    - `yaml_output` (used by `summary`): For generating YAML-formatted output.
- **External Libraries:**
    - No direct external library dependencies, although `fdf` might be considered an external utility for input parsing.
- **Interactions with other components:**
    - **`esl_system_m`:** `geometry_t` is a fundamental component of the main `system_t` type (typically as `system%geo`), making it central to the overall description of the simulation setup.
    - **Basis Set Generation (e.g., `esl_basis_m`):** The atomic coordinates, species types, and cell parameters from `geometry_t` are essential inputs for constructing the basis set for the calculation.
    - **Hamiltonian Construction & Solvers:** All parts of the code that build or operate on the Hamiltonian, calculate potentials, or solve for electronic states depend critically on the geometric information provided by `geometry_t`.
    - **Force Calculations (`esl_force_m`):** Atomic positions from `geometry_t` are needed to compute forces.
    - **Input Data Source:** The `init` routine acts as the primary interface for translating human-readable input (via FDF files) into the internal `geometry_t` representation.
    - **Output and Visualization:** The `summary` routine provides data that can be used for logging or as a basis for visualization tools.
