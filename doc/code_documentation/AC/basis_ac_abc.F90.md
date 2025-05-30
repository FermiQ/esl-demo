## Overview

The `esl_basis_ac_abc_t` module serves as a foundational module for defining atomic-centered (AC) orbital basis sets. It introduces three public derived types:
-   `orbital_ac_t`: Represents a single atomic orbital, detailing its quantum numbers (l, m), its radial part (via a pointer to a `pspiof_meshfunc_t` object), its effective cutoff radius, and its initial occupation.
-   `state_ac_t`: Groups `orbital_ac_t` objects that belong to a specific type of atomic site (typically a unique atomic species). It also stores the maximum cutoff radius among these orbitals.
-   `basis_ac_abc_t`: This is an `abstract` derived type that extends `basis_base_t` (from `esl_basis_base_m`). It's designed as an Abstract Base Class (ABC) for concrete AC basis set implementations. It aggregates all common data structures necessary for an AC basis, such as site information, mappings between sites and orbitals, the total number of orbitals, the sparsity pattern for matrices, and the overlap matrix.

The module's primary purpose is to provide a common structural definition for AC basis sets, potentially to manage complex dependencies or to enforce a standard interface for different AC basis implementations.

## Key Components

- **Module:** `esl_basis_ac_abc_t`
    - **Description:** Defines abstract base components and types for constructing atomic-centered (AC) basis sets.

### Public Derived Types

#### `orbital_ac_t`
- **Description:** Represents a single atomic-centered basis orbital. The full wavefunction $\psi(\mathbf{r}) = R(r) Y_l^m(\hat{\mathbf{r}})$ can be reconstructed from its components.
- **Components:**
    - `l (integer)`: The angular momentum quantum number of the orbital. Default: `0`.
    - `m (integer)`: The magnetic quantum number of the orbital. Default: `0`.
    - `r_cut (real(dp))`: The effective cutoff radius beyond which this orbital is considered negligible. Default: `0.0_dp`.
    - `occ (real(dp))`: The initial occupation of this orbital, typically derived from the configuration of the pseudo-atom used to generate the pseudopotential. Default: `0.0_dp`.
    - `R (type(pspiof_meshfunc_t), pointer :: R => null())`: A pointer to a `pspiof_meshfunc_t` object (from the `pspiof_m` library). This object stores the tabulated radial part of the basis function, $R(r)$.

#### `state_ac_t`
- **Description:** Represents a unique "state" or type of atomic site, typically corresponding to a specific atomic species (e.g., all Silicon atoms might be of the same `state_ac_t`). It groups all the basis orbitals that originate from this type of site.
- **Components:**
    - `n_orbital (integer)`: The total number of basis orbitals associated with this state type (e.g., for Si with s and p valence, this could be 1 for s + 3 for p = 4, if only one radial function per l is used; or more if multiple zetas/shapes per l are included). Default: `0`.
    - `r_cut (real(dp))`: The maximum effective cutoff radius (`r_cut`) among all `orbital_ac_t` objects belonging to this state. Default: `0.0_dp`.
    - `orb (type(orbital_ac_t), pointer :: orb(:) => null())`: A pointer to an array of `orbital_ac_t` objects that form the basis for this state type.

#### `basis_ac_abc_t`
- **Description:** An `abstract` derived type that serves as a base class for concrete atomic-centered basis set implementations. It extends `basis_base_t`, thereby inheriting components like `size` (total number of basis functions) and `grid` (an auxiliary real-space grid).
- **Components (in addition to those inherited from `basis_base_t`):**
    - `Q (real(dp))`: The total number of electrons, typically summed from the initial occupations (`occ`) of all basis orbitals. Default: `0.0_dp`.
    - `n_site (integer)`: The number of atomic sites (i.e., atoms) around which the basis functions are centered. Default: `0`.
    - `xyz (real(dp), allocatable :: xyz(:,:))`: A 2D allocatable array of shape `(3, n_site)` storing the Cartesian coordinates of each atomic site.
    - `n_state (integer)`: The number of unique site types (i.e., distinct `state_ac_t` definitions, usually one per atomic species). Default: `0`.
    - `state (type(state_ac_t), allocatable :: state(:))`: An allocatable array where each element is a `state_ac_t` object, describing the orbital composition for each unique site type.
    - `n_orbital (integer)`: The total number of atomic basis orbitals in the entire system. This is the sum of `state(is)%n_orbital` for each site, considering the species at that site. Default: `0`. (Note: this should be equal to the inherited `size`).
    - `site_state_idx (integer, allocatable :: site_state_idx(:))`: An array of length `n_site` that maps each site index to an index in the `state` array (i.e., identifies the species type of each atom).
    - `site_orbital_start (integer, allocatable :: site_orbital_start(:))`: An array of length `n_site + 1`. `site_orbital_start(i)` gives the global index of the first orbital centered on site `i`. `site_orbital_start(n_site+1)` would store `n_orbital + 1`.
    - `orbital_site (integer, allocatable :: orbital_site(:))`: An array of length `n_orbital` that maps each global orbital index back to the site index it is centered on.
    - `sparse_pattern (type(sparse_pattern_t))`: An object of type `sparse_pattern_t` (from `esl_sparse_pattern_m`) that defines the sparsity pattern (non-zero structure) for matrices like the overlap and Hamiltonian in this basis.
    - `S (type(sparse_matrix_t))`: An object of type `sparse_matrix_t` (from `esl_sparse_matrix_m`) used to store the overlap matrix $S_{ij} = \langle \psi_i | \psi_j \rangle$ in a sparse format.

## Important Variables/Constants
- The definitions of the three public types (`orbital_ac_t`, `state_ac_t`, `basis_ac_abc_t`) and their components are the main entities of this module. No public module-level parameters are defined.

## Usage Examples
```fortran
! TODO: Add usage example
! This module defines types, primarily an abstract base class.
! A concrete implementation would look like:
!
! module my_concrete_ac_basis_m
!   use esl_basis_ac_abc_t, only: basis_ac_abc_t
!   ! ... other imports ...
!   implicit none
!
!   type, extends(basis_ac_abc_t) :: my_ac_basis_t
!     ! ... potentially additional components specific to this implementation ...
!   contains
!     procedure, public :: init => my_init_ac_basis
!     ! ... other procedures implementing or overriding deferred procedures ...
!   end type my_ac_basis_t
!
! contains
!   subroutine my_init_ac_basis(this, ...)
!     class(my_ac_basis_t), intent(inout) :: this
!     ! ... arguments ...
!     ! Call parent's init if it had one, or initialize inherited components:
!     ! call this%basis_base_t%init(...)
!     ! this%Q = ...
!     ! allocate(this%xyz(...))
!     ! ... etc. ...
!   end subroutine my_init_ac_basis
!
! end module my_concrete_ac_basis_m
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec`: Provides precision kinds like `dp` (double precision).
    - `pspiof_m`: Provides the `pspiof_meshfunc_t` type, which is used in `orbital_ac_t` to store the radial part of basis functions. This links the basis set definition to pseudopotential data handling.
    - `esl_basis_base_m`: Provides `basis_base_t`, the parent type of `basis_ac_abc_t`, from which components like `size` (total number of orbitals) and `grid` (auxiliary real-space grid) are inherited.
    - `esl_sparse_pattern_m`: Provides `sparse_pattern_t`, used to define the non-zero structure of matrices.
    - `esl_sparse_matrix_m`: Provides `sparse_matrix_t`, used to store the overlap matrix `S`.
- **External Libraries:**
    - **PspIOf Library** (Indirectly via `pspiof_m`): The radial functions (`R` component of `orbital_ac_t`) are pointers to data structures managed by the PspIOf library, which reads them from pseudopotential files.
- **Interactions with other components:**
    - **Concrete Implementations (e.g., `esl_basis_ac_m`):** The `basis_ac_t` type defined in `esl_basis_ac_m` `extends` `basis_ac_abc_t`. This means `basis_ac_t` inherits all the data components declared in `basis_ac_abc_t` and then provides concrete implementations for procedures that manipulate this data (like `init`, `get_psi`, etc.).
    - **`esl_species_m`:** Data originating from `species_t` objects (such as the radial parts of pseudo-atomic orbitals, their l-quantum numbers, and occupations) are used by the `init` routine of a concrete AC basis implementation to populate the `orbital_ac_t` and `state_ac_t` structures.
    - **Hamiltonian and Density Matrix Modules for AC basis (e.g., `esl_hamiltonian_ac_m`, `esl_density_matrix_ac_m`):** These modules operate on instances of a concrete AC basis type. They utilize the components defined in `basis_ac_abc_t` (like `n_orbital`, `sparse_pattern`, `S`, and individual orbital data via `state(:)%orb(:)`) for their calculations.
    - **Modularity and Extensibility:** The use of an abstract base class like `basis_ac_abc_t` can facilitate code organization, allowing for different types of atomic-centered basis sets to be developed in the future while adhering to a common set of data structures and (potentially) a common procedural interface. The comment about "circular dependencies" suggests this structure helps in managing the relationships between different modules defining and using the basis set.
</tbody>
</table>
