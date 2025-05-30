## Overview

The `esl_hamiltonian_pw_m` module is responsible for managing and operating with the electronic Hamiltonian within a plane-wave (PW) basis set. It defines a derived type `hamiltonian_pw_t` that primarily holds a pointer to a `potential_t` object (which contains the Hartree, exchange-correlation, and external local potentials).

The core functionality of this module revolves around its `eigensolver` procedure, which uses an iterative method, specifically the ELSI RCI (Rayleigh Quotient Iteration or similar, from the `elsi_rci` module), to find the lowest eigenvalues and corresponding eigenvectors (wavefunctions) of the Kohn-Sham Hamiltonian. This iterative solver repeatedly applies the Hamiltonian to trial wavefunctions.

The module provides routines for:
-   Applying the full Hamiltonian ($H\psi = [-\frac{1}{2}\nabla^2 + V_{local}]\psi$) to a wavefunction (`hamiltonian_pw_apply`). The kinetic term is applied in reciprocal space, while the local potential term is applied in real space using Fast Fourier Transforms (FFTs).
-   Applying only the local part of the potential ($V_{local}\psi$) (`hamiltonian_pw_apply_local`).
-   Applying a simple preconditioner (`hamiltonian_preconditioner`), often used to accelerate convergence in iterative eigensolvers.

An internal helper type `work_matrix_t` is defined to facilitate the management of temporary matrices required by the ELSI RCI solver.

## Key Components

- **Module:** `esl_hamiltonian_pw_m`
    - **Description:** Manages Hamiltonian operations and eigenvalue solutions for plane-wave basis sets, primarily interfacing with the ELSI RCI solver.

- **Type (Internal):** `work_matrix_t`
    - **Description:** A simple container for a 2D allocatable complex matrix, used for workspace arrays in the `eigensolver`.
    - **Components:**
        - `mat (complex(dp), allocatable :: mat(:, :))`: The complex matrix.

- **Type:** `hamiltonian_pw_t`
    - **Description:** Represents the plane-wave Hamiltonian, primarily through its link to the local potential.
    - **Components:**
        - `pot (type(potential_t), pointer :: pot => null())`: A pointer to a `potential_t` object (from `esl_potential_m`), which provides the Hartree, exchange-correlation, and external local potential components defined on a real-space grid.
    - **Procedures (Public, Type-Bound):**
        - `init(this, pot, comm)`: Initializes the Hamiltonian object.
        - `eigensolver(this, pw, states)`: Solves the Kohn-Sham equations to find eigenvalues and eigenvectors.
    - **Finalizer (Type-Bound):**
        - `finalizer(this)`: Nullifies the `pot` pointer.

- **Public Module Subroutines:**
    - `hamiltonian_pw_apply(pot_in, pw_basis_in, psi_in, hpsi_out)`: Applies the full Hamiltonian ($T + V_{local}$) to a wavefunction `psi_in`, result in `hpsi_out`.
    - `hamiltonian_pw_apply_local(pw_basis_in, pot_in, psi_in, hpsi_out)`: Applies only the local potential part ($V_{local}$) to `psi_in`, adding to `hpsi_out`.
    - `hamiltonian_preconditioner(pw_basis_in, psi_in, hpsi_out)`: Applies a diagonal preconditioner in reciprocal space to `psi_in`, result in `hpsi_out`.

### Selected Procedure Details

#### `init(this, pot_in, comm)` (Type-Bound)
- **Description:** Initializes the `hamiltonian_pw_t` object `this`. Its primary action is to associate the internal pointer `this%pot` with the provided `pot_in` object (of type `potential_t`). The `comm` (MPI communicator) argument is accepted but not directly used in the current body, though it might have been intended for ELSI initialization steps that are now handled elsewhere.
- **Arguments:**
    - `this (class(hamiltonian_pw_t), intent(inout))`.
    - `pot_in (type(potential_t), target, intent(in))`: The potential object containing $V_H, V_{XC}, V_{ext}$.
    - `comm (integer, intent(in))`: MPI communicator.

#### `eigensolver(this, pw, states)` (Type-Bound)
- **Description:** This is the core routine for finding eigenvalues and eigenvectors. It employs an iterative solver strategy driven by `rci_omm` from the `elsi_rci` module.
    1.  Sets up workspace: Allocates an array `work` of `work_matrix_t` (size 28) to hold various matrices required by `rci_omm`. Copies initial guess wavefunctions from `states%states(:,:,1)%zcoef` into `work(1)%mat`.
    2.  **Iterative Loop:** Enters a `do` loop that calls `rci_omm`. `rci_omm` returns a `task` identifier and information (`iS`) about which matrices in `work` to operate on.
    3.  A `select case (task)` block performs the requested operation:
        *   `ELSI_RCI_CONVERGE`: Loop terminates.
        *   `ELSI_RCI_H_MULTI`: Computes $B = H \cdot A$. Calls `hamiltonian_pw_apply` for columns of $A$ (from `work(iS%Aidx)%mat`) storing results in $B$ (`work(iS%Bidx)%mat`).
        *   `ELSI_RCI_S_MULTI`: Computes $B = S \cdot A$. For an orthonormal PW basis ($S=I$), this is $B=A$.
        *   `ELSI_RCI_P_MULTI`: Computes $B = P \cdot A$ (preconditioner application). Calls `hamiltonian_preconditioner`.
        *   Other tasks (`ELSI_RCI_GEMM`, `_AXPY`, `_COPY`, `_TRACE`, `_DOT`, `_SCALE`) perform corresponding BLAS/LAPACK-like matrix operations on specified `work` matrices.
    4.  **Final Diagonalization (Rayleigh-Ritz):** After convergence, `rci_omm` will have typically produced an orthonormal basis $A$ (e.g., `work(1)%mat`) that spans the desired subspace and the projected Hamiltonian $H_{proj} = A^H H A$ (e.g., in `work(21)%mat`). This projected Hamiltonian is diagonalized using `zheev` (LAPACK).
    5.  The eigenvalues from `zheev` are stored in `states%eigenvalues`.
    6.  The final eigenvectors are reconstructed as $C_{final} = A \cdot C_{subspace}$ (where $C_{subspace}$ are eigenvectors from `zheev`) via `zgemm` and stored in `states%states(:,:,1)%zcoef`.
- **Arguments:**
    - `this (class(hamiltonian_pw_t), intent(in))`.
    - `pw (type(basis_pw_t), intent(in))`: The plane-wave basis definition.
    - `states (type(states_t), intent(inout))`: Electronic states object; updated with new eigenvalues and eigenvectors.

#### `hamiltonian_pw_apply(pot_in, pw_basis_in, psi_in, hpsi_out)` (Module Subroutine)
- **Description:** Applies the full Kohn-Sham Hamiltonian ($H = -\frac{1}{2}\nabla^2 + V_{local}$) to an input wavefunction `psi_in` (PW coefficients), producing `hpsi_out`.
    1.  **Kinetic Energy ($T\psi$):** In reciprocal space, this is $\frac{1}{2}|\mathbf{G}+\mathbf{k}|^2 \psi(\mathbf{G}+\mathbf{k})$. This is applied element-wise: `hpsi_out(ic) = -0.5_dp * psi_in(ic) * pw_basis_in%gmod2(ic)`. (Note: `gmod2` usually stores $|\mathbf{G}+\mathbf{k}|^2$, so the sign might be for $H\psi = E\psi$, or it's simply $0.5 |\mathbf{G}+\mathbf{k}|^2$).
    2.  **Local Potential ($V_{local}\psi$):** Calls `hamiltonian_pw_apply_local` to calculate the local potential part and add it to `hpsi_out`.
- **Arguments:**
    - `pot_in (type(potential_t), intent(in))`: Provides the local potentials.
    - `pw_basis_in (type(basis_pw_t), intent(in))`: Provides PW basis info like `gmod2`.
    - `psi_in (complex(dp), intent(in) :: psi_in(:))`: Input wavefunction coefficients.
    - `hpsi_out (complex(dp), intent(inout) :: hpsi_out(:))`: Output $(H\psi)$ coefficients.

#### `hamiltonian_preconditioner(pw_basis_in, psi_in, hpsi_out)` (Module Subroutine)
- **Description:** Applies a simple diagonal preconditioner in reciprocal space. For each G-vector component `ic`: `hpsi_out(ic) = 2.0_dp * psi_in(ic) / (pw_basis_in%gmod2(ic) + 1.0e-08_dp)`. This is related to $(T + \text{const})^{-1}$.
- **Arguments:** Similar to `hamiltonian_pw_apply` but `pot_in` is not needed.

#### `hamiltonian_pw_apply_local(pw_basis_in, pot_in, psi_in, hpsi_out)` (Module Subroutine)
- **Description:** Applies the local potential part of the Hamiltonian, $V_{local}\psi$, where $V_{local} = V_{Hartree} + V_{XC} (+ V_{External}$, often included in $V_{Hartree}$ from Poisson solver).
    1.  Transforms `psi_in` (PW coefficients) to its real-space representation `psi_rs` using `pw2grid` (from `esl_utils_pw_m`).
    2.  In real space, multiplies `psi_rs(ip)` by `(pot_in%hartree(ip) + pot_in%vxc(ip))` for each grid point `ip`, storing the result in `vpsi_rs(ip)`.
    3.  Transforms `vpsi_rs` back to reciprocal space using `grid2pw` (from `esl_utils_pw_m`) and adds the result to `hpsi_out` (using `add=.true.` in `grid2pw`).
- **Arguments:** Similar to `hamiltonian_pw_apply`.

#### `finalizer(this)` (Type-Bound)
- **Description:** Nullifies the `this%pot` pointer.
- **Arguments:** `this (type(hamiltonian_pw_t), intent(inout))`.

## Important Variables/Constants
- **`pot` (pointer to `potential_t`)**: The key component of `hamiltonian_pw_t`, providing access to local potential terms.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage within an SCF module:
!
! type(hamiltonian_pw_t) :: H_pw
! type(potential_t), target :: current_potential
! type(basis_pw_t) :: pw_basis
! type(states_t) :: current_states
! integer :: mpi_communicator
!
! ! ... (Initialize current_potential, pw_basis, current_states, mpi_communicator) ...
!
! call H_pw%init(current_potential, mpi_communicator)
!
! ! In SCF loop:
! ! ... (current_potential is updated based on new density) ...
! call H_pw%eigensolver(pw_basis, current_states) ! Updates current_states
!
! ! H_pw is finalized automatically when it goes out of scope.
```

## Dependencies and Interactions

- **Core ESL Modules:**
    - `prec`: For precision kinds (`dp`, `ip`).
    - `esl_basis_pw_m`: Provides `basis_pw_t` (G-vector info, grid for FFTs).
    - `esl_grid_m`: (Indirectly via `basis_pw_t`) Defines the real-space grid.
    - `esl_potential_m`: Provides `potential_t` (Hartree, XC, external potentials).
    - `esl_states_m`: Provides `states_t` for initial guess wavefunctions and storing results.
    - `esl_utils_pw_m`: Provides `pw2grid` and `grid2pw` for FFT-based transformations between real and reciprocal space.
- **External Libraries/Modules (Crucial):**
    - **ELSI RCI (via `elsi_rci`, `elsi_rci_constants`, `elsi_rci_omm`, `elsi_rci_precision` modules):** The `eigensolver` is fundamentally dependent on the `rci_omm` routine from this library to perform the iterative diagonalization.
    - **BLAS/LAPACK:** Routines like `zgemm`, `zaxpy` (BLAS Level 1 and 3) and `zheev` (LAPACK) are used by the `eigensolver` for various matrix manipulations and the final subspace diagonalization.
- **Interactions with other components:**
    - **Main Hamiltonian Wrapper (`esl_hamiltonian_m`):** `hamiltonian_pw_t` is the specialized type for PW calculations and would be a component of the main `hamiltonian_t` (e.g., `main_H%pw`).
    - **SCF Cycle (`esl_scf_m`):** The `eigensolver` method of `hamiltonian_pw_t` is the workhorse called in each SCF iteration after the electronic potential (managed by `this%pot`) has been updated based on the previous density.
    - **Potential Module (`esl_potential_m`):** The `hamiltonian_pw_t` holds a pointer to a `potential_t` object, from which it gets the local potential components ($V_H, V_{XC}, V_{ext}$) needed for the $V_{local}\psi$ operation.
    - **Basis Module (`esl_basis_pw_m`):** Provides essential plane-wave basis information like G-vector magnitudes (`gmod2`), mappings (`gmap`), and the associated real-space grid for FFTs.
    - **Utility Modules (`esl_utils_pw_m`):** These provide the FFT wrappers (`pw2grid`, `grid2pw`) used in applying the local potential.
</tbody>
</table>
