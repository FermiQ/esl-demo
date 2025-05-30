## Overview

The `esl_basis_ac_m` module defines the `basis_ac_t` derived type, which represents an atomic-centered (AC) orbital basis set for electronic structure calculations. This type extends `esl_basis_ac_abc_t`, an abstract base class for atomic-centered bases, inheriting its core structure and likely some procedures.

The primary responsibilities of this module include:
-   Initializing the AC basis set by processing atomic geometry and species information. This involves reading pseudo-atomic orbitals from pseudopotential files (via `esl_species_m` which uses `pspiof_m`) and constructing basis functions as products of these radial parts and spherical harmonics (Y_lm).
-   Managing the properties of these basis functions, such as their angular momentum (l, m), radial parts, cutoff radii, and initial occupations (derived from the pseudo-atom configuration).
-   Organizing orbitals by atomic sites and species types.
-   Generating a sparsity pattern (`sparse_pattern`) appropriate for matrices (like overlap and Hamiltonian) in this AC basis, based on the spatial extent of the orbitals.
-   Calculating the overlap matrix (`S`) using this sparsity pattern.
-   Providing methods to evaluate the basis functions ($\psi_i(\mathbf{r})$) at any given point in space (`get_psi`).
-   Offering a utility to construct an initial density matrix based on atomic orbital occupations (`atomic_density_matrix`).
-   It also maintains an auxiliary real-space grid, initialized based on a specified mesh cutoff energy.

## Key Components

- **Module:** `esl_basis_ac_m`
    - **Description:** Manages atomic-centered orbital basis sets, including their construction from pseudopotential data, overlap calculation, and evaluation.
- **Type:** `basis_ac_t`, extends `esl_basis_ac_abc_t`
    - **Description:** Represents an atomic-centered orbital basis set. It inherits components from `basis_ac_abc_t`, such as:
        - `grid (type(grid_t))`: An auxiliary real-space grid.
        - `n_site (integer)`: Number of atomic sites (atoms).
        - `xyz(3, n_site) (real(dp), allocatable)`: Cartesian coordinates of these sites.
        - `site_state_idx(n_site) (integer, allocatable)`: Maps each site to a "state" index (representing a unique species type).
        - `n_state (integer)`: Number of unique species types.
        - `state(n_state) (type(state_ac_t), allocatable)`: An array where each element (`state_ac_t`) describes the set of basis orbitals originating from one unique atomic species. `state_ac_t` likely contains:
            - `n_orbital (integer)`: Number of basis orbitals for this species type.
            - `r_cut (real(dp))`: Maximum cutoff radius among orbitals of this species type.
            - `orb(n_orbital) (type(orbital_ac_t), allocatable)`: Array of `orbital_ac_t`.
        - `orbital_ac_t` (likely nested or imported type) would contain:
            - `l (integer)`: Angular momentum quantum number.
            - `m (integer)`: Magnetic quantum number.
            - `R (type(pspiof_meshfunc_t), pointer)`: Pointer to the tabulated radial part of the basis function (from PspIOf).
            - `occ (real(dp))`: Initial occupation of this orbital from the pseudo-atom.
            - `r_cut (real(dp))`: Effective cutoff radius for this specific orbital.
        - `n_orbital (integer)`: Total number of atomic orbitals in the basis set across all sites.
        - `site_orbital_start(n_site+1) (integer, allocatable)`: For each site, the starting index in the global list of orbitals. `site_orbital_start(n_site+1)` stores `n_orbital + 1`.
        - `orbital_site(n_orbital) (integer, allocatable)`: Maps each global orbital index back to its atomic site index.
        - `Q (real(dp))`: Total electronic charge from the sum of initial atomic orbital occupations.
        - `size (integer)`: Total number of basis functions (same as `n_orbital`).
        - `sparse_pattern (type(sparse_pattern_t))`: Defines the sparsity structure for matrices based on orbital overlaps.
        - `S (type(sparse_matrix_t))`: The overlap matrix $S_{ij} = \langle \psi_i | \psi_j \rangle$.
    - **Type-Bound Procedures (Public):**
        - `init(this, geo)`: Initializes the basis set.
        - `summary(this)`: Prints a detailed summary.
        - `get_psi` (generic interface for `get_psi_all`, `get_psi_single`): Evaluates basis functions.
        - `atomic_density_matrix(this, DM)`: Constructs a density matrix from atomic occupations.
    - **Finalizer (Type-Bound):**
        - `finalizer(this)`: Deallocates all allocatable components, with special care for shared radial function pointers.

### Selected Type-Bound Procedures

#### `init(this, geo)`
- **Description:** Initializes the `basis_ac_t` object using geometry information (`geo` of type `geometry_t`).
    1.  Reads `MeshCutOff` (for auxiliary grid) and `Basis.AC.FR.Cutoff` (for radial function extent) from FDF input.
    2.  Initializes `this%grid` using `ndims_from_meshcutoff` (a private module subroutine that determines grid dimensions based on `MeshCutOff` and cell vectors).
    3.  Copies atomic site coordinates and species indices from `geo` into `this%xyz` and `this%site_state_idx`.
    4.  Iterates over each unique atomic species defined in `geo%species`:
        *   Determines the number of basis orbitals derived from this species by summing `2*l+1` for each pseudo-atomic radial orbital provided in the species' pseudopotential data.
        *   Allocates `this%state(is)%orb` to hold these orbitals.
        *   For each pseudo-atomic radial orbital (e.g., 2s, 2p, 3d from the pseudopotential file of species `is`):
            *   Retrieves its angular momentum `l`, its tabulated radial function (`mesh_R` of `pspiof_meshfunc_t`), and its occupation `occ` in the pseudo-atom (using `geo%species(is)%get_radial_orbital`).
            *   Determines its effective cutoff radius `r_cut` (using `geo%species(is)%get_radial_orbital_rmax` with `Basis.AC.FR.Cutoff`).
            *   Updates `this%state(is)%r_cut` to be the maximum of such cutoffs for the current species.
            *   For each `m` from `-l` to `l`, it creates an `orbital_ac_t` entry: sets `l`, `m`, associates its radial part `R` as a pointer to `mesh_R`, sets occupation per m-channel, and stores `r_cut`. The same `mesh_R` pointer is used for all m-channels of a given (l, radial_function) pair to save memory.
    5.  Calculates the total number of orbitals `this%n_orbital` and the total initial electronic charge `this%Q` by summing occupations.
    6.  Populates mapping arrays: `this%orbital_site` (orbital index to site index) and `this%site_orbital_start` (site index to first orbital index on that site). Sets `this%size = this%n_orbital`.
    7.  Calls `create_sparse_pattern_ac_create(this, this%sparse_pattern)` to determine the sparsity pattern of matrices based on orbital overlaps (using `r_cut` values).
    8.  Calls `overlap_matrix_ac_calculate(this, this%sparse_pattern, this%S)` to compute the overlap matrix `S`.
- **Arguments:**
    - `this (class(basis_ac_t), intent(inout))`: The AC basis object to initialize.
    - `geo (type(geometry_t), intent(in))`: The system geometry, providing atomic species, their pseudopotential data, and coordinates.

#### `atomic_density_matrix(this, DM)`
- **Description:** Constructs a sparse density matrix `DM` representing the initial atomic guess (sum of isolated atom densities). The diagonal elements corresponding to each basis orbital are filled with the initial occupation of that orbital (from `this%state(is)%orb(iio)%occ`). If the input `DM` is not already initialized with a sparse pattern, a minimal diagonal sparse pattern is created.
- **Arguments:**
    - `this (class(basis_ac_t), intent(in))`: The AC basis object.
    - `DM (type(sparse_matrix_t), intent(inout))`: The sparse density matrix to be populated.

#### `get_psi` (Generic Interface for `get_psi_all` and `get_psi_single`)
- **`get_psi_all(this, state_idx, r_rel, psi_values)`**: Calculates all basis functions `psi_values(:)` associated with a given `state_idx` (species type index) at a point `r_rel` (vector from orbital center to evaluation point). `psi_values` must be large enough.
- **`get_psi_single(this, state_idx, orbital_in_state_idx, r_rel, psi_value)`**: Calculates a single basis function `psi_value` corresponding to the `orbital_in_state_idx`-th orbital of species type `state_idx`, at point `r_rel`.
- **Mechanism:** Both routines use `grylmr` (from `esl_numeric_m`) to get the Y_lm angular part and `pspiof_meshfunc_eval` (from `pspiof_m`) to evaluate the tabulated radial part (`orb%R`) at the given distance. The final basis function value is their product.

#### `finalizer(this)`
- **Description:** Deallocates all allocatable arrays within `this`. Special care is taken for the radial function pointers (`this%state(is)%orb(io)%R`), as multiple m-orbitals might share the same radial function data. The logic attempts to identify unique radial function pointers and deallocate each only once before nullifying all pointers.
- **Arguments:** `this (type(basis_ac_t))`: The AC basis object to finalize.

## Important Variables/Constants
- All components of `basis_ac_t` are crucial for defining the atomic-centered basis, especially the `state` array (which describes orbitals per species type) and the resulting `n_orbital`, `sparse_pattern`, and overlap matrix `S`.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage:
!
! use esl_basis_ac_m
! use esl_geometry_m
! type(basis_ac_t) :: ac_basis
! type(geometry_t) :: system_geo
!
! ! ... (Initialize system_geo, e.g., by reading from input) ...
!
! ! Initialize the atomic-centered basis
! call ac_basis%init(system_geo)
!
! ! Print a summary of the basis
! call ac_basis%summary()
!
! ! ... (ac_basis can now be used by Hamiltonian and density matrix modules) ...
!
! ! ac_basis will be finalized automatically when it goes out of scope.
```

## Dependencies and Interactions

- **Base Type:** `esl_basis_ac_abc_t`: `basis_ac_t` extends this abstract base class.
- **Core ESL Modules:**
    - `prec`: For precision kinds (`dp`).
    - `fdf`: For reading configuration parameters (e.g., `MeshCutOff`, `Basis.AC.FR.Cutoff`).
    - `yaml_output`: Used by `summary` for formatted output.
    - `esl_numeric_m`: Provides `grylmr` for calculating spherical harmonics (used in `get_psi`).
    - `esl_message_m`: For error reporting.
    - `esl_species_m`: Provides `species_t` and its methods (like `get_radial_orbital`, `get_radial_orbital_rmax`) which are fundamental for extracting radial parts of basis functions from pseudopotential data.
    - `esl_geometry_m`: Provides `geometry_t`, which is the primary input for `init`, supplying atomic coordinates and species information.
    - `esl_grid_m`: The `grid_t` type (from `esl_basis_ac_abc_t`'s `grid` component) is initialized and used.
- **AC Specific ESL Modules:**
    - `esl_create_sparse_pattern_ac_m`: Provides `create_sparse_pattern_ac_create` used during `init` to determine the sparsity structure of matrices.
    - `esl_overlap_matrix_ac_m`: Provides `overlap_matrix_ac_calculate` used during `init` to compute the overlap matrix `S`.
    - `esl_sparse_matrix_m`: Provides `sparse_matrix_t` (for `S` and for `DM` in `atomic_density_matrix`).
- **External Libraries:**
    - `pspiof_m` (PspIOf): This library interface is used via `esl_species_m` to read and interpret pseudopotential files, from which the radial parts of the basis functions are derived. `pspiof_meshfunc_t` objects store these radial parts.
- **Interactions with other components:**
    - **System Setup (`esl_system_m`, `esl_main_m`):** If an atomic-centered basis is selected for a calculation, an instance of `basis_ac_t` would be created and initialized, typically as part of the overall `system_t` or `basis_t` object.
    - **Hamiltonian Construction (`esl_hamiltonian_ac_m`):** This module would heavily rely on `basis_ac_t` to:
        - Obtain the overlap matrix `S`.
        - Get the sparsity pattern for the Hamiltonian matrix.
        - Evaluate matrix elements which involve integrals of basis functions with parts of the Hamiltonian operator. This often requires evaluating basis functions (using `get_psi`).
    - **Density Matrix Calculation (`esl_density_ac_m`, `esl_density_matrix_ac_m`):** These modules use the `basis_ac_t` to understand the structure of the density matrix, its size, and its sparsity. The `atomic_density_matrix` method provides an initial guess.
    - **Input Configuration (`fdf`):** Cutoff parameters influencing the basis set construction and auxiliary grid are read from FDF input files.
</tbody>
</table>
