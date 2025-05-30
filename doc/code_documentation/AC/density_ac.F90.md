## Overview

The `esl_density_ac_m` module provides the `density_ac_t` derived type, specifically designed for managing electron density representations when using atomic-centered (AC) orbital basis sets. This type extends `density_base_t` (which contributes the real-space electron density `rho` on a grid) and adds components for handling the density matrix (`DM`) and the energy-density matrix (`EDM`) in the sparse AC basis.

The module offers a suite of functionalities including:
-   Initialization of the density object, including allocating space for `rho` and initializing the sparse `DM` and `EDM` based on a sparsity pattern provided by the AC basis.
-   Generation of an initial guess for both the density matrix (based on atomic valence electron configurations) and the corresponding real-space density.
-   Calculation of the real-space density `rho` from the current density matrix `DM`.
-   Calculation of a new density matrix `DM` and energy-density matrix `EDM` from a given set of electronic states (wavefunction coefficients, occupations, and eigenvalues).
-   A static module procedure `add_density_matrix` to project a given sparse density matrix onto the real-space grid density `rho`.
-   Calculation and printing of Mulliken population analysis via `mulliken_summary`.
-   Computation of the residual (difference) between two density matrices, often used for checking SCF convergence.

## Key Components

- **Module:** `esl_density_ac_m`
    - **Description:** Manages electron density for atomic-centered basis sets, including real-space density and sparse density matrices.

- **Type:** `density_ac_t`, extends `density_base_t`
    - **Description:** A derived type for AC density representations.
    - **Inherited Components (from `density_base_t`):**
        - `rho (real(dp), allocatable :: rho(:))`: Stores the electron density on a real-space grid.
    - **Specific Components:**
        - `DM (type(sparse_matrix_t))`: The density matrix $P_{\mu\nu}$ in the atomic orbital basis, stored as a sparse matrix.
        - `EDM (type(sparse_matrix_t))`: The energy-density matrix $E_{\mu\nu}$ (often $P_{\mu\nu} \times \epsilon_\nu$ or similar), also sparse.
    - **Type-Bound Procedures (Public):**
        - `init(this, basis)`: Initializes the density object.
        - `residue(this, other)`: Calculates the max absolute difference between `this%DM` and `other%DM`.
        - `guess(this, basis)`: Generates an initial guess for `DM` and `rho`.
        - `calculate(this, basis)`: Computes `this%rho` from `this%DM`.
        - `calculate_density_matrix(this, states)`: Computes `this%DM` and `this%EDM` from a `states_t` object.
    - **Finalizer (Type-Bound):**
        - `finalizer(this)`: Deallocates `rho` and finalizes `DM` and `EDM`.

### Module Subroutines (Public)
- **`add_density_matrix(basis, DM_in, rho_out)`**: (Static procedure) Adds the contribution of a given sparse density matrix `DM_in` to a real-space density array `rho_out`, using basis functions from `basis`. $\rho(\mathbf{r}) = \sum_{\mu,\nu} \psi_\mu(\mathbf{r}) DM_{\mu\nu} \psi_\nu(\mathbf{r})$.
- **`mulliken_summary(this_density_ac, basis, S_matrix)`**: Calculates and prints Mulliken charges based on `this_density_ac%DM` and the overlap matrix `S_matrix`.

### Selected Procedure Details

#### `init(this, basis)`
- **Description:** Initializes a `density_ac_t` object. It allocates the real-space density array `this%rho` based on the number of points in `basis%grid%np`. It also initializes the sparse matrix objects `this%DM` and `this%EDM` using the sparsity pattern provided by `basis%sparse_pattern`.
- **Arguments:**
    - `this (class(density_ac_t), intent(inout))`.
    - `basis (type(basis_ac_t), intent(in))`: The AC basis set definition.

#### `guess(this, basis)`
- **Description:** Generates an initial guess for the density.
    1.  Calls `basis%atomic_density_matrix(this%DM)` to populate `this%DM` with values corresponding to a superposition of neutral atom electronic configurations.
    2.  Initializes `this%EDM%M` (the non-zero elements of EDM) to zeros.
    3.  Calls the module procedure `add_density_matrix(basis, this%DM, this%rho)` to compute the real-space density `this%rho` from this initial atomic `DM`. `this%rho` is zeroed before the call.
- **Arguments:**
    - `this (class(density_ac_t), intent(inout))`.
    - `basis (type(basis_ac_t), intent(in))`.

#### `calculate(this, basis)`
- **Description:** Calculates the real-space density `this%rho` from the current density matrix `this%DM`. It first zeros out `this%rho` and then calls `add_density_matrix(basis, this%DM, this%rho)` to project the density matrix onto the grid.
- **Arguments:**
    - `this (class(density_ac_t), intent(inout))`.
    - `basis (type(basis_ac_t), intent(in))`.

#### `calculate_density_matrix(this, states)`
- **Description:** Calculates the density matrix `this%DM` and energy-density matrix `this%EDM` from a `states_t` object (which contains wavefunction coefficients, occupations, and eigenvalues). The formulas are (schematically, for orbital indices $\mu, \nu$ and state index $kbs$ summing over k-points, spins, and bands):
    -   $DM_{\mu\nu} = \sum_{kbs} C_{\mu,kbs}^* \text{occ}_{kbs} C_{\nu,kbs}$
    -   $EDM_{\mu\nu} = \sum_{kbs} C_{\mu,kbs}^* \text{occ}_{kbs} \epsilon_{kbs} C_{\nu,kbs}$
    (where $C$ are wavefunction coefficients, $\text{occ}$ are occupation numbers, $\epsilon$ are eigenvalues, and k-point weights are included in $\text{occ}_{kbs}$). The current implementation handles real coefficients; complex coefficients have a "to be implemented" note.
- **Arguments:**
    - `this (class(density_ac_t), intent(inout))`.
    - `states (type(states_t), intent(in))`: The electronic states information.

#### `residue(this, other) result(res)`
- **Description:** Computes the residual between two density matrices, defined as the maximum absolute difference between corresponding elements of `this%DM%M` and `other%DM%M`. This is often used to check for SCF convergence.
- **Arguments:**
    - `this (class(density_ac_t), intent(in))`.
    - `other (class(density_ac_t), intent(in))`.
- **Returns:** `res (real(dp))`: The maximum absolute difference.

#### `finalizer(this)`
- **Description:** Type-bound final procedure. Deallocates `this%rho` if allocated, and calls `this%DM%delete()` and `this%EDM%delete()` to finalize the sparse matrix objects.
- **Arguments:** `this (type(density_ac_t), intent(inout))`.

## Important Variables/Constants
- **`DM (type(sparse_matrix_t))`**: The sparse density matrix in the AC basis.
- **`EDM (type(sparse_matrix_t))`**: The sparse energy-density matrix.
- **`rho (real(dp), allocatable :: rho(:))`**: (Inherited) The electron density on a real-space grid.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage within an SCF cycle:
!
! type(density_ac_t) :: density_obj
! type(basis_ac_t) :: ac_basis
! type(states_t) :: current_electron_states
! ! ... (initialize ac_basis, current_electron_states) ...
!
! call density_obj%init(ac_basis)
! call density_obj%guess(ac_basis) ! Initial guess for DM and rho
!
! ! In SCF loop:
! ! ... (current_electron_states are updated by eigensolver) ...
! call density_obj%calculate_density_matrix(current_electron_states) ! Update DM and EDM from states
! call density_obj%calculate(ac_basis) ! Update rho from the new DM
! ! ... (use density_obj%rho for next potential calculation, check convergence with residue) ...
```

## Dependencies and Interactions

- **Base Type:** `esl_density_base_m`: `density_ac_t` extends `density_base_t`.
- **Core ESL Modules:**
    - `prec`: For precision kinds (`dp`).
    - `esl_basis_ac_m`: Provides `basis_ac_t`, which is crucial for defining the AC basis, evaluating basis functions (used in `add_density_matrix`), and providing the sparsity pattern for `DM` and `EDM`.
    - `esl_grid_m`: The `grid_t` (accessed via `basis%grid`) defines the real-space mesh for `rho`.
    - `esl_sparse_matrix_m`, `esl_sparse_pattern_m`: Provide types and methods for sparse matrices (`DM`, `EDM`) and their patterns.
    - `esl_states_m`: Provides `states_t`, used as input for `calculate_density_matrix`.
    - `yaml_output` (used by `mulliken_summary`).
- **Interactions with other components:**
    - **`esl_density_m` (Main Density Wrapper):** `density_ac_t` is the specialized type used by `esl_density_m` when calculations are performed with an AC basis. The `density_t` object in `esl_density_m` would contain an instance of `density_ac_t`.
    - **SCF Cycle (e.g., `esl_scf_m`):**
        - `guess` provides the initial density.
        - `calculate_density_matrix` is called after new electronic states are found to update `DM` and `EDM`.
        - `calculate` is then called to get the updated `rho` from the new `DM`.
        - `residue` is used to compare density matrices for SCF convergence.
    - **Potential Calculation (e.g., `esl_potential_m`):** The real-space density `this%rho` is a primary input for calculating Hartree and exchange-correlation potentials.
    - **Mulliken Analysis (`mulliken_summary`):** This routine uses the `DM` and an overlap matrix `S` (from `basis_ac_t` or elsewhere) to compute Mulliken charges.
    - **ELSI Interface (Indirectly):** While this module itself might not directly call ELSI for DM/EDM calculation (that's often done by `esl_calc_density_matrix_ac_m` which then populates `DM`/`EDM`, or the states are solved by ELSI), the `states_t` object it consumes can be a result of ELSI calculations.
</tbody>
</table>
