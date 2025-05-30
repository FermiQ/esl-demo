## Overview

The `esl_hamiltonian_ac_m` module is responsible for constructing, managing, and diagonalizing the Hamiltonian matrix within the context of an atomic-centered (AC) orbital basis set. It defines the `hamiltonian_ac_t` derived type, which stores the different components of the Hamiltonian—kinetic energy, non-local pseudopotential (Kleinman-Bylander projectors), and the total Self-Consistent Field (SCF) Hamiltonian—as sparse matrices.

The module provides a comprehensive suite of procedures for:
-   Initializing the sparse matrix structures for these Hamiltonian components.
-   Calculating the matrix elements for the non-SCF dependent parts:
    -   Kinetic energy ($-\frac{1}{2}\nabla^2$) via `hamiltonian_ac_laplacian`.
    -   Non-local Kleinman-Bylander pseudopotential contributions via `hamiltonian_ac_Vkb` (a private module routine called by `calculate_H0`).
-   Assembling the initial SCF Hamiltonian ($H_0$) from these fixed components.
-   Adding matrix elements of a local potential (e.g., sum of external, Hartree, and XC potentials) to the SCF Hamiltonian, handled by `hamiltonian_ac_potential`.
-   Solving the generalized eigenvalue problem $(H - \epsilon S)C = 0$ to obtain eigenvalues and eigenvectors (electronic states). The current implementation includes a solver for the Gamma point using real arithmetic (`dsygv` from LAPACK).

The calculations of matrix elements leverage the grid and basis function evaluation capabilities provided by `esl_grid_m` and `esl_basis_ac_m`.

## Key Components

- **Module:** `esl_hamiltonian_ac_m`
    - **Description:** Manages the construction and diagonalization of the Hamiltonian for atomic-centered basis sets using sparse matrix representations.

- **Public Module Subroutines (Matrix Element Calculators):**
    - `hamiltonian_ac_laplacian(basis, H_kin_out)`: Calculates the kinetic energy matrix elements $\langle\psi_i | -\frac{1}{2}\nabla^2 | \psi_j\rangle$ and stores them in the sparse matrix `H_kin_out`.
    - `hamiltonian_ac_potential(basis, V_on_grid, H_V_out)`: Calculates the matrix elements $\langle\psi_i | V(\mathbf{r}) | \psi_j\rangle$ for a given local potential `V_on_grid` (defined on the real-space grid) and adds them to the sparse matrix `H_V_out`.

- **Type:** `hamiltonian_ac_t`
    - **Description:** A derived type to store and manage different sparse matrix components of the Hamiltonian in an AC basis.
    - **Components:**
        - `kin (type(sparse_matrix_t))`: Sparse matrix storing the kinetic energy term ($T_{\mu\nu}$).
        - `vkb (type(sparse_matrix_t))`: Sparse matrix storing the contribution from non-local Kleinman-Bylander pseudopotential projectors.
        - `H (type(sparse_matrix_t), allocatable :: H(:))`: An allocatable array of sparse matrices, representing the total SCF Hamiltonian for each spin component (e.g., `H(1)` for spin-up or non-spin-polarized).
    - **Procedures (Public, Type-Bound):**
        - `init(this, sparse_pattern, nspin)`: Initializes the `kin`, `vkb`, and `H` sparse matrices.
        - `calculate_H0(this, basis, geom)`: Computes the non-SCF dependent terms (`kin` and `vkb`).
        - `setup_H0(this)`: Assembles the initial Hamiltonian $H = T + V_{KB}$.
        - `add_potential(this, basis, V_on_grid)`: Adds the matrix elements of a local potential `V_on_grid` to `this%H(1)`.
        - `eigensolver(this, basis, states)`: Solves the generalized eigenvalue problem.
    - **Finalizer (Type-Bound):**
        - `finalizer(this)`: Deallocates the `kin`, `vkb`, and `H` sparse matrices.

### Selected Procedure Details

#### `init(this, sparse_pattern, nspin)` (Type-Bound)
- **Description:** Initializes the `hamiltonian_ac_t` object. It initializes `this%kin` and `this%vkb` as sparse matrices using the provided `sparse_pattern`. It then allocates `this%H` for `nspin` components and initializes each `this%H(ispin)` also using `sparse_pattern`.
- **Arguments:**
    - `this (class(hamiltonian_ac_t), intent(inout))`.
    - `sparse_pattern (type(sparse_pattern_t), intent(in), target)`: The sparsity pattern for all Hamiltonian components.
    - `nspin (integer, intent(in))`: Number of spin components.

#### `calculate_H0(this, basis, geom)` (Type-Bound)
- **Description:** Calculates the non-SCF (geometry-dependent but density-independent) parts of the Hamiltonian.
    1.  Zeros out `this%kin%M` and calls `hamiltonian_ac_laplacian(basis, this%kin)` to fill the kinetic energy matrix.
    2.  Zeros out `this%vkb%M` and calls `hamiltonian_ac_Vkb(geom, basis, this%vkb)` (a private module subroutine) to fill the non-local pseudopotential matrix.
- **Arguments:**
    - `this (class(hamiltonian_ac_t), intent(inout))`.
    - `basis (type(basis_ac_t), intent(in))`: AC basis definition.
    - `geom (type(geometry_t), intent(in))`: System geometry.

#### `setup_H0(this)` (Type-Bound)
- **Description:** Constructs the initial Hamiltonian(s) before the SCF cycle by summing the pre-calculated kinetic and non-local pseudopotential terms: `this%H(ispin)%M = this%kin%M + this%vkb%M` for each spin component.
- **Arguments:** `this (class(hamiltonian_ac_t), intent(inout))`.

#### `add_potential(this, basis, V_on_grid)` (Type-Bound)
- **Description:** Adds the matrix elements of a given local potential `V_on_grid` (which is defined on the real-space grid) to the current SCF Hamiltonian. It calls the module subroutine `hamiltonian_ac_potential(basis, V_on_grid, this%H(1))`. Currently, it adds only to the first spin component `this%H(1)`.
- **Arguments:**
    - `this (class(hamiltonian_ac_t), intent(inout))`.
    - `basis (type(basis_ac_t), intent(in))`: AC basis definition.
    - `V_on_grid (real(dp), intent(in) :: V(:))`: The local potential values on the real-space grid.

#### `eigensolver(this, basis, states)` (Type-Bound)
- **Description:** Solves the generalized eigenvalue problem $H C = \epsilon S C$ for the current SCF Hamiltonian `this%H` and overlap `basis%S`.
    - If `states%complex_states` is true, it calls an internal `eig_k()` routine (currently a stub that prints "to be implemented").
    - If `states%complex_states` is false (real arithmetic, Gamma point), it calls an internal `eig_gamma()` routine:
        - For each spin component:
            - It converts the sparse `this%H(ispin)` and `basis%S` matrices to dense matrices (`H_dense`, `S_dense`).
            - It calls the LAPACK routine `dsygv` (option `1`, `'V'`, `'U'`) to solve for all eigenvalues and eigenvectors.
            - Checks `info` from `dsygv` for errors.
            - Copies the lowest `states%nstates` eigenvalues into `states%eigenvalues(:,ispin,1)`.
            - Copies the corresponding eigenvectors (columns of `H_dense` after `dsygv`) into `states%states(:,ispin,1)%dcoef(:)`.
- **Arguments:**
    - `this (class(hamiltonian_ac_t), intent(in))`.
    - `basis (type(basis_ac_t), intent(in))`: AC basis, providing the overlap matrix `S`.
    - `states (type(states_t), intent(inout))`: The electronic states object where eigenvalues and eigenvectors (coefficients) will be stored.

#### `finalizer(this)` (Type-Bound)
- **Description:** Deallocates the sparse matrix objects `this%kin`, `this%vkb`, and all components of `this%H(:)` by calling their respective `delete()` methods.
- **Arguments:** `this (type(hamiltonian_ac_t), intent(inout))`.

## Important Variables/Constants
- **`kin`, `vkb`, `H`**: Instances of `sparse_matrix_t` that store the key Hamiltonian components.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage within an SCF cycle:
!
! type(hamiltonian_ac_t) :: H_ac
! type(basis_ac_t) :: ac_basis
! type(geometry_t) :: sys_geometry
! type(states_t) :: current_states
! type(sparse_pattern_t) :: sp_pattern
! integer :: n_spin_comp
! real(dp), allocatable :: local_potential_on_grid(:)
!
! ! ... (Initialize ac_basis, sys_geometry, current_states, sp_pattern, n_spin_comp) ...
! ! ... (Allocate and fill local_potential_on_grid) ...
!
! call H_ac%init(ac_basis%sparse_pattern, n_spin_comp) ! ac_basis%sparse_pattern should be used
! call H_ac%calculate_H0(ac_basis, sys_geometry)      ! Calculate T and V_NL
!
! ! Inside SCF loop:
! call H_ac%setup_H0()                                 ! H = T + V_NL
! call H_ac%add_potential(ac_basis, local_potential_on_grid) ! H = H + V_local
! call H_ac%eigensolver(ac_basis, current_states)      ! Solve for new states
!
! ! H_ac is finalized automatically when it goes out of scope.
```

## Dependencies and Interactions

- **Core ESL Modules:**
    - `prec`: For precision kinds (`dp`).
    - `fdf`: For reading parameters like `Basis.AC.KB.Cutoff`.
    - `esl_basis_ac_m`: Provides `basis_ac_t` (basis functions, overlap matrix `S`, grid access, orbital info).
    - `esl_geometry_m`: Provides `geometry_t` (atom positions, species info for projectors).
    - `esl_states_m`: Provides `states_t` for storing results of `eigensolver`.
    - `esl_grid_m`: The `grid_t` (from `basis%grid`) is used extensively for evaluating basis functions, projectors, and integrating matrix elements.
    - `esl_sparse_matrix_m`, `esl_sparse_pattern_m`: For `sparse_matrix_t` and `sparse_pattern_t` types and methods.
    - `esl_constants_m`: Used in `hamiltonian_ac_laplacian` for `PI`.
- **External Libraries/Modules:**
    - `pspiof_m`: Used indirectly via `basis` and `geom` to access pseudopotential projector information (`geom%species(:)%get_projector`, `geom%species(:)%get_projector_rmax`) and radial parts of basis functions.
    - **LAPACK:** The `dsygv` routine is used in `eigensolver` for solving the dense generalized eigenvalue problem.
- **Interactions with other components:**
    - **Main Hamiltonian Wrapper (`esl_hamiltonian_m`):** `hamiltonian_ac_t` is the specialized type for AC calculations, likely a component of the main `hamiltonian_t` (e.g., `main_H%ac`).
    - **SCF Cycle (`esl_scf_m`):** This module is central to the SCF loop for AC basis sets.
        - `init` and `calculate_H0` are called during setup.
        - In each iteration: `setup_H0` re-initializes H from fixed parts, `add_potential` adds the current density-dependent local potential (from `esl_potential_m`), and then `eigensolver` finds the new electronic states.
    - **Potential Module (`esl_potential_m`):** The local potential (sum of external, Hartree, and XC) computed by `potential_t` is passed as `V_on_grid` to `this%add_potential`.
    - **Basis Module (`esl_basis_ac_m`):** Provides all necessary information about the basis functions, their representation on the grid, the overlap matrix `S`, and the sparsity pattern.
</tbody>
</table>
