## Overview

The `esl_create_sparse_pattern_ac_m` module provides a specialized subroutine, `create_sparse_pattern_ac_create`, designed to determine and construct a sparsity pattern for matrices that arise in calculations using atomic-centered (AC) basis sets. The sparsity pattern indicates which matrix elements are potentially non-zero, allowing for efficient storage and computation using sparse matrix techniques.

The fundamental principle employed is that a matrix element between two basis orbitals, $\psi_i$ and $\psi_j$, is considered potentially non-zero only if the orbitals exhibit spatial overlap. This module approximates this condition by examining the distance between the atomic centers on which the orbitals are located and comparing it to the sum of their effective cutoff radii. If the distance is within this sum, the corresponding matrix element $(i,j)$ (and its symmetric counterpart $(j,i)$) is marked as a non-zero entry in the sparsity pattern. Diagonal elements and off-diagonal elements between orbitals on the same atom are also included.

## Key Components

- **Module:** `esl_create_sparse_pattern_ac_m`
    - **Description:** Generates a sparsity pattern for matrices based on the spatial extent (cutoff radii) and positions of atomic orbitals in an atomic-centered (AC) basis.

### Public Subroutine

#### `create_sparse_pattern_ac_create(basis, sp)`
- **Description:** This subroutine populates a `sparse_pattern_t` object (`sp`) based on the geometric arrangement and cutoff radii of atomic orbitals defined in the `basis` (of type `basis_ac_abc_t`).
    1.  **Reset Pattern:** It begins by calling `sp%delete()` to clear any existing data in the `sp` object.
    2.  **Initialization:** It determines `max_no` (the maximum number of orbitals on any single atom) and `no` (the total number of orbitals in the `basis`). The `sp` object is then initialized using `sp%init(no, no, np=max_no * 20)`, where `np` is an initial estimate for the number of non-zero elements per row, used for pre-allocation (the heuristic `max_no * 20` is used).
    3.  **Pattern Generation Loop:** It iterates through all pairs of atomic sites (`ia`, `ja`) where basis functions are centered:
        *   **Intra-atomic elements (ia == ja):** It first calls an internal helper subroutine `add_elements(ia, ia, 0.0_dp)`. This helper (detailed below) adds entries to `sp` for all diagonal elements $(io, io)$ and all unique off-diagonal elements $(io, jo)$ where both orbitals `io` and `jo` are centered on the same atom `ia`.
        *   **Inter-atomic elements (ja > ia):** To process each pair of distinct sites once, it loops for `ja` from `ia + 1` to `basis%n_site`.
            *   It retrieves a general cutoff radius for the species at site `ia` (`ir_cut_species = basis%state(is)%r_cut`) and site `ja` (`jr_cut_species = basis%state(js)%r_cut`).
            *   It calculates the distance `dist` between the centers of atoms `ia` and `ja`.
            *   If `dist <= ir_cut_species + jr_cut_species` (a broad check using species-level maximum cutoffs), it calls the internal helper `add_elements(ia, ja, dist)`. This helper then performs more detailed checks using individual orbital cutoff radii.
    4.  **Finalization:** After iterating through all relevant site pairs, `call sp%finalize()` is invoked. This step typically sorts the collected non-zero indices, removes any duplicates, and builds the final compressed representation (e.g., CSR or CSC pointers) of the sparsity pattern within the `sp` object.

- **Internal Helper Subroutine:** `add_elements(ia, ja, dist)`
    -   **Case 1: Intra-atomic (`ia == ja`)**
        -   Iterates through all orbitals `io` centered on atom `ia`.
        -   Adds the diagonal element `(io, io)` to `sp`.
        -   Iterates through orbitals `jo` (where `jo > io`) also on atom `ia` and adds `(io, jo)` and `(jo, io)` to `sp`.
    -   **Case 2: Inter-atomic (`ia < ja`, ensured by caller context or internal check)**
        -   Iterates through each orbital `io` on site `ia` and each orbital `jo` on site `ja`.
        -   Retrieves the specific cutoff radius for orbital `io` (`ir_cut_orb = basis%state(is)%orb(local_io_idx)%r_cut`) and orbital `jo` (`jr_cut_orb`).
        -   If the distance `dist` (between atomic centers `ia` and `ja`) is less than or equal to the sum of these specific orbital cutoff radii (`ir_cut_orb + jr_cut_orb`), it adds entries `(io, jo)` and `(jo, io)` to the sparsity pattern `sp`. This is a more refined check than the species-level cutoff check done by the caller.

- **Arguments:**
    - `basis (class(basis_ac_abc_t), intent(in))`: The atomic-centered basis set definition. This provides access to atomic site coordinates (`basis%xyz`), the species type at each site (`basis%site_state_idx`), and for each species, the list of its orbitals along with their individual cutoff radii (`basis%state(:)%orb(:)%r_cut`). It also uses the species-level maximum cutoff (`basis%state(:)%r_cut`).
    - `sp (type(sparse_pattern_t), intent(inout))`: The sparsity pattern object that will be populated by this routine.
- **Returns:** Not applicable (Subroutine).

## Important Variables/Constants
- This module does not define any public Fortran `parameter` constants. The logic relies on the cutoff radii defined within the `basis` object.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage (typically called from a basis set initialization routine):
!
! module esl_basis_ac_m
!   use esl_basis_ac_abc_t
!   use esl_create_sparse_pattern_ac_m, only: create_sparse_pattern_ac_create
!   use esl_sparse_pattern_m, only: sparse_pattern_t
!   ! ... other imports ...
!
!   type, extends(basis_ac_abc_t) :: basis_ac_t
!     ! ... components inherited and potentially new ones ...
!   contains
!     procedure, public :: init => init_basis_ac
!     ! ... other procedures ...
!   end type basis_ac_t
!
! contains
!   subroutine init_basis_ac(this, geo)
!     class(basis_ac_t), intent(inout) :: this
!     type(geometry_t), intent(in) :: geo ! From esl_geometry_m
!
!     ! ... (initialize orbital definitions, site info, xyz, r_cut in 'this' from 'geo') ...
!
!     ! Create the sparsity pattern for matrices in this basis
!     ! 'this%sparse_pattern' is of type(sparse_pattern_t), inherited from basis_ac_abc_t
!     call create_sparse_pattern_ac_create(this, this%sparse_pattern)
!
!     ! ... (now this%sparse_pattern can be used to initialize sparse matrices like overlap S) ...
!   end subroutine init_basis_ac
!
! end module esl_basis_ac_m
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only: dp`: Imports the `dp` kind specifier for double precision real numbers.
    - `esl_basis_ac_abc_t, only: basis_ac_abc_t`: This is crucial as the input `basis` argument is of this type (or a type that extends it). The routine accesses various components of `basis_ac_abc_t` such as site coordinates, species information, orbital lists, and cutoff radii.
    - `esl_sparse_pattern_m, only: sparse_pattern_t`: Provides the `sparse_pattern_t` derived type and its associated methods (`delete`, `init`, `add`, `finalize`) which are used to construct the sparsity pattern.
- **External Libraries:**
    - This module does not have direct dependencies on external libraries.
- **Interactions with other components:**
    - **`esl_basis_ac_m` (Concrete AC Basis Implementation):** The `init` routine of `basis_ac_t` (defined in `esl_basis_ac_m`) is a primary caller of `create_sparse_pattern_ac_create`. After setting up the orbital definitions, it calls this routine to determine the `sparse_pattern` member of the `basis_ac_t` object.
    - **Sparse Matrix Construction (e.g., `esl_overlap_matrix_ac_m`, `esl_hamiltonian_ac_m`):** Modules responsible for building matrices (like the overlap matrix or Hamiltonian matrix) in the AC basis will use the `sparse_pattern_t` object generated by this module. The pattern dictates which matrix elements need to be computed and stored, allowing for efficient initialization and use of `sparse_matrix_t` objects.
    - **Linear Algebra Operations:** Downstream, sparse linear algebra libraries or routines (e.g., eigensolvers like ELSI, or iterative solvers) will rely on this sparsity pattern for efficient operations, as they only need to consider the non-zero elements of the matrices.
    - **Performance:** The accuracy and tightness of the cutoff radii used to determine sparsity can significantly impact both the memory footprint of sparse matrices and the computational cost of matrix operations. Overly generous cutoffs lead to denser patterns, while overly tight cutoffs could erroneously zero out small but significant matrix elements.
</tbody>
</table>
