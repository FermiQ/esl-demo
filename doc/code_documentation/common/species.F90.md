## Overview

The `esl_species_m` module defines a derived type `species_t`, which serves as a container for all information related to a specific atomic species within an electronic structure calculation. A primary role of this module is to interface with the `pspiof_m` library (PspIOf) to read and parse data from pseudopotential files. Pseudopotentials are used to simplify calculations by replacing the strong core electron potential and core electrons with a weaker effective potential, allowing treatment of only valence electrons.

A `species_t` object stores the species label (e.g., "Si", "O"), the complete pseudopotential data structure from PspIOf, the valence charge density, the local part of the pseudopotential, information about non-local projectors, details of pseudo-atomic radial orbitals, the ionic charge (charge of the core), and the number of valence electrons. The module provides methods to initialize a species by reading its pseudopotential file and to access its various components, including functions to determine effective cutoff radii for orbitals and projectors.

## Key Components

- **Module:** `esl_species_m`
    - **Description:** Manages data for atomic species, primarily by reading and interpreting pseudopotential files using the PspIOf library.
- **Type:** `species_t`
    - **Description:** A derived type holding all data associated with a single atomic species, largely derived from its pseudopotential file.
    - **Components:**
        - `label (character(len=10))`: A short, descriptive label for the species (e.g., "H", "Si").
        - `psp (type(pspiof_pspdata_t))`: An object from the PspIOf library that stores the entire dataset read from the pseudopotential file.
        - `rho (type(pspiof_meshfunc_t))`: A PspIOf mesh function object representing the tabulated valence charge density of the pseudo-atom.
        - `vlocal (type(pspiof_meshfunc_t))`: A PspIOf mesh function object representing the tabulated local part of the pseudopotential.
        - `n_projectors (integer)`: The number of non-local projectors defined in the pseudopotential.
        - `projectors (type(pspiof_projector_t), allocatable, private :: projectors(:))`: A private allocatable array storing the individual non-local projector function objects from PspIOf.
        - `n_radial_orbitals (integer)`: The number of pseudo-atomic radial orbitals (valence states) defined.
        - `radial_orbitals (type(pspiof_state_t), allocatable, private :: radial_orbitals(:))`: A private allocatable array storing the pseudo-atomic radial orbital objects (including their wavefunctions and occupations) from PspIOf.
        - `z_ion (real(dp))`: The ionic charge (i.e., charge of the nucleus + core electrons; Z_core).
        - `q (real(dp))`: The number of valence electrons for this species.
    - **Procedures (Public, Type-Bound):**
        - `init`: Initializes the species by reading and parsing a pseudopotential file.
        - `get_radial_orbital`: Retrieves details (l-quantum number, radial wavefunction, occupation) of a specific pseudo-atomic radial orbital.
        - `get_radial_orbital_rmax`: Calculates an effective maximum radius for a specified radial orbital.
        - `get_projector`: Retrieves details (l-quantum number, energy, radial part) of a specific non-local projector.
        - `get_projector_rmax`: Calculates an effective maximum radius for a specified projector.
        - `summary`: Prints a brief summary of the species information.
    - **Finalizer (Type-Bound):**
        - `finalizer`: Frees resources allocated by the PspIOf library for this species.

### Type-Bound Procedures

#### `init(this, label, filename)`
- **Description:** Initializes a `species_t` object by reading data from a specified pseudopotential file using the PspIOf library.
    1.  Sets the species `this%label`.
    2.  Allocates `this%psp` (a `pspiof_pspdata_t` structure).
    3.  Reads the pseudopotential data from the file specified by `filename` into `this%psp`. Reports an error if the file cannot be read.
    4.  Extracts and stores key information from `this%psp`:
        *   Valence charge density into `this%rho`.
        *   Pseudo-atomic radial orbitals into `this%radial_orbitals` (allocating the array based on `pspiof_pspdata_get_n_states`).
        *   Ionic charge (`this%z_ion`) and number of valence electrons (`this%q`).
        *   The local potential `this%vlocal`. This involves a workaround: it retrieves the mesh and a `pspiof_potential_t` object, evaluates the potential on the mesh points, and then initializes `this%vlocal` (a `pspiof_meshfunc_t`) with this tabulated data.
        *   Non-local projectors into `this%projectors` (allocating based on `pspiof_pspdata_get_n_projectors`).
- **Arguments:**
    - `this (class(species_t), intent(inout))`: The `species_t` object to be initialized.
    - `label (character(len=*), intent(in))`: The label for this atomic species (e.g., "Si").
    - `filename (character(len=*), intent(in))`: The path to the pseudopotential file.

#### `get_radial_orbital(this, io, ll, radial_orbital_mf, occ)`
- **Description:** Retrieves information about the `io`-th pseudo-atomic radial orbital stored in the species data.
- **Arguments:**
    - `this (class(species_t), intent(in))`: The species object.
    - `io (integer, intent(in))`: The 1-based index of the desired radial orbital.
    - `ll (integer, optional, intent(out))`: If present, assigned the angular momentum quantum number (l) of the orbital.
    - `radial_orbital_mf (type(pspiof_meshfunc_t), optional, intent(out))`: If present, assigned the `pspiof_meshfunc_t` object representing the radial part of the orbital's wavefunction.
    - `occ (real(dp), optional, intent(out))`: If present, assigned the occupation number of this orbital in the reference atomic configuration.

#### `get_radial_orbital_rmax(this, io, tolerance) result(rmax)`
- **Description:** Calculates an effective maximum radius (`rmax`) for the `io`-th pseudo-atomic radial orbital. This radius is determined by finding the point on the radial mesh beyond which the absolute value of the orbital's wavefunction multiplied by the radius (`abs(psi(r)*r)`) drops below a given `tolerance`.
- **Arguments:**
    - `this (class(species_t), intent(in))`: The species object.
    - `io (integer, intent(in))`: The 1-based index of the radial orbital.
    - `tolerance (real(dp), intent(in))`: The cutoff tolerance value.
- **Returns:** `rmax (real(dp))`: The calculated effective maximum radius.

#### `get_projector(this, ip, projector_mf, energy, l)`
- **Description:** Retrieves information about the `ip`-th non-local projector function. Similar to the local potential extraction in `init`, this routine reconstructs a `pspiof_meshfunc_t` object for the projector by evaluating the `pspiof_projector_t` on the radial mesh, as a direct getter might not be available in PspIOf for this conversion.
- **Arguments:**
    - `this (class(species_t), intent(in))`: The species object.
    - `ip (integer, intent(in))`: The 1-based index of the desired projector.
    - `projector_mf (type(pspiof_meshfunc_t), intent(out))`: The `pspiof_meshfunc_t` object representing the radial part of the projector.
    - `energy (real(dp), intent(out))`: The energy associated with this projector (often related to the atomic state it's derived from).
    - `l (integer, intent(out))`: The angular momentum quantum number (l) of the projector.

#### `get_projector_rmax(this, ip, tolerance) result(rmax)`
- **Description:** Calculates an effective maximum radius (`rmax`) for the `ip`-th non-local projector function. This radius is determined by finding the point on the radial mesh beyond which the absolute value of the projector function itself drops below a given `tolerance`.
- **Arguments:**
    - `this (class(species_t), intent(in))`: The species object.
    - `ip (integer, intent(in))`: The 1-based index of the projector.
    - `tolerance (real(dp), intent(in))`: The cutoff tolerance value.
- **Returns:** `rmax (real(dp))`: The calculated effective maximum radius.

#### `summary(this)`
- **Description:** Prints a brief summary of the species to standard output using YAML format. The summary includes the species `label`, its ionic charge (`z_ion`), and its electronic valence charge (`q`).
- **Arguments:**
    - `this (class(species_t), intent(in))`: The species object to summarize.

#### `finalizer(this)`
- **Description:** This is a type-bound `final` subroutine automatically invoked when a `species_t` object is deallocated. It calls `pspiof_meshfunc_free(this%vlocal)` to free the workaround-generated local potential mesh function and `pspiof_pspdata_free(this%psp)` to release all resources allocated by the PspIOf library for the main pseudopotential data.
- **Arguments:**
    - `this (type(species_t))`: The `species_t` object being finalized.

## Important Variables/Constants
- **`label`**: User-friendly species identifier.
- **`psp`**: Holds all raw pseudopotential data via PspIOf.
- **`rho`, `vlocal`**: Tabulated valence density and local potential.
- **`projectors`, `radial_orbitals`**: Arrays holding non-local projectors and pseudo-atomic orbitals.
- **`z_ion`, `q`**: Core ionic charge and number of valence electrons, crucial for charge neutrality and setup.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage:
!
! use esl_species_m
! type(species_t) :: silicon_species
! character(len=80) :: psp_file_path = "Si.psf" ! Example pseudopotential file
! integer :: l_orb
! type(pspiof_meshfunc_t) :: orbital_func
!
! ! Initialize the silicon species by reading its pseudopotential file
! call silicon_species%init("Si", psp_file_path)
!
! ! Print a summary of the species
! call silicon_species%summary()
!
! ! Get information about the first radial orbital
! if (silicon_species%n_radial_orbitals > 0) then
!   call silicon_species%get_radial_orbital(1, ll=l_orb, radial_orbital=orbital_func)
!   print *, "First radial orbital has l = ", l_orb
!   ! ... can now use orbital_func which is a pspiof_meshfunc_t ...
!   call pspiof_meshfunc_free(orbital_func) ! Important to free if it was allocated by get_radial_orbital
! end if
!
! ! silicon_species will be finalized automatically when it goes out of scope.
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `esl_grid_m`: While not directly `use`d for a `grid_t` object within `species_t` itself, the underlying `pspiof_m` library extensively uses mesh/grid concepts (`pspiof_mesh_t`, `pspiof_meshfunc_t`) which are conceptually similar to grids managed by `esl_grid_m`.
    - `esl_message_m`: Used for error reporting via `message_error` if pseudopotential file reading fails.
    - `prec`: Provides precision kinds like `dp`.
    - `pspiof_m`: This is the **core external library interface** used by `esl_species_m`. All operations related to reading, parsing, and extracting data from pseudopotential files are delegated to functions from `pspiof_m` (e.g., `pspiof_pspdata_alloc`, `pspiof_pspdata_read`, `pspiof_pspdata_get_rho_valence`, etc.).
    - `yaml_output`: Used by the `summary` procedure for YAML-formatted output.
- **External Libraries:**
    - **PspIOf Library:** The library accessed via the `pspiof_m` Fortran interface module. This library is essential for handling pseudopotential files.
- **Interactions with other components:**
    - **`esl_geometry_m`:** The `geometry_t` type (from `esl_geometry_m`) contains an array of `species_t` objects (`geometry%species(:)`). The `geometry_t%init` routine is responsible for creating and initializing these `species_t` objects based on the species defined in the main input file.
    - **`esl_potential_m`:** The `potential_t%compute_ext_loc` method (in `esl_potential_m`) uses the `this%vlocal` component (the local pseudopotential) of each `species_t` object to construct the total external local potential for the system.
    - **Hamiltonian Construction (e.g., `esl_hamiltonian_ac_m`, `esl_hamiltonian_pw_m`):** These modules would use the `get_projector` method to retrieve the non-local projector functions required for building the non-local part of the Kohn-Sham Hamiltonian. The `get_radial_orbital` method might be used for constructing atomic-orbital basis sets or for generating initial guesses for wavefunctions.
    - **Basis Set Generation:** In codes employing atomic orbitals as basis functions, the `get_radial_orbital` method could be directly used to obtain the radial parts of these basis functions from the pseudopotential definition.
    - **Input System:** The `init` routine relies on a filename, implying that the choice of pseudopotential for each species is part of the user's input to the overall application.
</tbody>
</table>
