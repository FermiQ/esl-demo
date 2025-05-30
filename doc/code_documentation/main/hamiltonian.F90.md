## Overview

The `esl_hamiltonian_m` module provides a high-level interface for managing and solving the electronic Hamiltonian of a system. It acts as a dispatcher, choosing between Plane-Wave (PW) specific routines (from `esl_hamiltonian_pw_m`) or Atomic-Centered (AC) orbital specific routines (from `esl_hamiltonian_ac_m`) based on the system's configured basis set. Its main functions are to initialize the Hamiltonian components, including the electronic potential, and to invoke the appropriate eigensolver to determine the electronic states (eigenvalues and eigenvectors).

## Key Components

- **Module:** `esl_hamiltonian_m`
    - **Description:** Serves as a wrapper for basis-dependent Hamiltonian operations, delegating tasks to either plane-wave or atomic-centered Hamiltonian modules.
- **Type:** `hamiltonian_t`
    - **Description:** A derived type that aggregates all data structures related to the Hamiltonian. This includes the electronic potential (`potential_t`) and specific types for plane-wave (`hamiltonian_pw_t`) and atomic-centered (`hamiltonian_ac_t`) Hamiltonians.
    - **Components:**
        - `potential (type(potential_t))`: Stores the total electronic potential experienced by the electrons.
        - `pw (type(hamiltonian_pw_t))`: Contains data and procedures specific to plane-wave Hamiltonian calculations.
        - `ac (type(hamiltonian_ac_t))`: Contains data and procedures specific to atomic-centered orbital Hamiltonian calculations.
    - **Procedures (public):**
        - `init`: Initializes the Hamiltonian, including its potential component and any basis-dependent fixed terms.
        - `eigensolver`: Solves the Kohn-Sham eigenvalue problem for the current Hamiltonian.

### Subroutine/Function: `init(this, basis, geo, states, periodic)`

- **Description:** Initializes the `hamiltonian_t` object. This is a crucial setup step that performs the following:
    1.  Initializes the `this%potential` object (of type `potential_t`) using the provided basis, geometry, states (for information like number of spins), and periodicity.
    2.  Based on `basis%type`:
        *   If `PLANEWAVES`: Calls `this%pw%init`, providing it with the newly initialized potential and the MPI communicator (`MPI_COMM_WORLD` if `WITH_MPI` is defined, otherwise a dummy value).
        *   If `ATOMCENTERED`:
            *   Calls `this%ac%init`, passing the sparse matrix pattern (obtained from `basis%ac%sparse_pattern`) and the number of electron spins (`states%nspin`). This sets up the structure for the AC Hamiltonian matrix.
            *   Calls `this%ac%calculate_H0(basis%ac, geo)` to compute and store the parts of the Hamiltonian that do not change during the Self-Consistent Field (SCF) cycle. These typically include the kinetic energy operator and the non-local pseudopotential projectors (Kleinman-Bylander projectors), which depend only on the atomic positions and the basis set.
- **Arguments:**
    - `this (class(hamiltonian_t))`: The Hamiltonian object to be initialized.
    - `basis (type(basis_t), intent(in))`: The basis set configuration (determines PW or AC path).
    - `geo (type(geometry_t), intent(in))`: The system's atomic geometry.
    - `states (type(states_t), intent(in))`: Information about electronic states (e.g., number of spins).
    - `periodic (logical, intent(in))`: A flag indicating whether the system has periodic boundary conditions.
- **Returns:** Not applicable (Subroutine).

### Subroutine/Function: `eigensolver(this, basis, states)`

- **Description:** This subroutine is responsible for solving the eigenvalue problem (e.g., Kohn-Sham equations) for the current Hamiltonian to obtain the electronic eigenvalues and eigenvectors (wavefunctions). It delegates the task to the appropriate basis-specific eigensolver:
    - If `basis%type` is `PLANEWAVES`: Calls `this%pw%eigensolver(basis%pw, states)`.
    - If `basis%type` is `ATOMCENTERED`: Calls `this%ac%eigensolver(basis%ac, states)`.
The `states` object is updated in place with the results.
- **Arguments:**
    - `this (class(hamiltonian_t), intent(inout))`: The Hamiltonian object. The potential component (`this%potential`) is assumed to be up-to-date before this call.
    - `basis (type(basis_t), intent(in))`: The basis set configuration.
    - `states (type(states_t), intent(inout))`: The electronic states object. On input, it may provide initial guesses or information like occupation numbers. On output, it is updated with the calculated eigenvalues and eigenvectors.
- **Returns:** Not applicable (Subroutine).

## Important Variables/Constants

- This module does not define its own public parameters but relies on `PLANEWAVES` and `ATOMCENTERED` from the `esl_basis_m` module (accessed via `basis%type`) to direct its workflow.

## Usage Examples

```fortran
! TODO: Add usage example
! Assuming 'ham' is a hamiltonian_t object,
! 'sys_basis' is an initialized basis_t object,
! 'sys_geo' is an initialized geometry_t object,
! 'sys_states' is an initialized states_t object,
! 'is_periodic' is a logical flag.

! Initialize Hamiltonian
! call ham%init(sys_basis, sys_geo, sys_states, is_periodic)

! Update potential part of ham (ham%potential) based on new density (not shown here)
! ...

! Solve for eigenstates
! call ham%eigensolver(sys_basis, sys_states)
! The 'sys_states' object now contains the new eigenvalues and eigenvectors.
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only : dp, ip`: Imports double precision (`dp`) and integer precision (`ip`) type kinds.
    - `esl_basis_m`: Provides `basis_t` (including `basis%pw`, `basis%ac`, `basis%ac%sparse_pattern`) and the basis type identifiers `PLANEWAVES` and `ATOMCENTERED`.
    - `esl_geometry_m`: Provides `geometry_t` (atomic coordinates and cell information), used during initialization.
    - `esl_grid_m`: [How it's used - Not directly invoked by `esl_hamiltonian_m` itself, but likely used by its components, e.g., `esl_potential_m` for real-space grids or `esl_hamiltonian_pw_m` for k-space grids.]
    - `esl_hamiltonian_pw_m`: Provides `hamiltonian_pw_t` and its methods (`init`, `eigensolver`) for plane-wave based calculations.
    - `esl_hamiltonian_ac_m`: Provides `hamiltonian_ac_t` and its methods (`init`, `calculate_H0`, `eigensolver`) for atomic-centered orbital based calculations.
    - `esl_sparse_pattern_m`: Provides `sparse_pattern_t`, used by `esl_hamiltonian_ac_m` via `basis%ac%sparse_pattern` for efficient storage of sparse AC Hamiltonian matrices.
    - `esl_potential_m`: Provides `potential_t`. An instance of `potential_t` is a member of `hamiltonian_t` and is initialized by it. The potential is fundamental to defining the Hamiltonian.
    - `esl_states_m`: Provides `states_t`, which is used to pass information (like number of spins) to `init` and to receive the results (eigenvalues, eigenvectors) from `eigensolver`.
    - `mpi` (conditional, if `WITH_MPI` is defined): Used for MPI parallelization, specifically `MPI_COMM_WORLD` is passed to `esl_hamiltonian_pw_m%init`.
- **External Libraries:**
    - `MPI`: If compiled `WITH_MPI`, this library is used for distributed parallel computation, particularly relevant for the plane-wave parts.
- **Interactions with other components:**
    - **SCF Cycle:** This module is a critical part of the Self-Consistent Field (SCF) procedure.
        - `init` is called once at the beginning (or after geometry changes) to set up fixed parts of H.
        - In each SCF iteration:
            - The `potential` component of `hamiltonian_t` is updated using the latest electron density (from `esl_density_m`). This update is typically handled by procedures within `esl_potential_m` called by the main SCF driver.
            - `eigensolver` is then called to solve the Kohn-Sham equations with this updated potential.
    - **`esl_potential_m`:** The `hamiltonian_t` type contains and manages a `potential_t` object. The accuracy of this potential (derived from the electron density) directly affects the Hamiltonian.
    - **`esl_density_m`:** The electron density computed by `esl_density_m` is used to construct/update the potential in `this%potential`.
    - **`esl_states_m`:** The `eigensolver` method populates the `states_t` object with the computed eigenvalues and eigenvectors, which are then used, for example, to calculate the next iteration's electron density.
    - **Basis-specific modules (`esl_hamiltonian_pw_m`, `esl_hamiltonian_ac_m`):** `esl_hamiltonian_m` acts as a facade, delegating most of the complex work to these specialized modules.
