## Overview

The `esl_scf_m` module is central to the electronic structure calculation, as it implements the Self-Consistent Field (SCF) procedure. The SCF method iteratively refines the electronic density and potential of a system until a self-consistent solution is achieved, yielding the ground state electronic structure. This module manages the entire SCF loop, including the initialization of relevant quantities (densities, Hamiltonian), the iterative solution of the effective single-particle equations (e.g., Kohn-Sham equations), application of density mixing techniques to aid convergence, and checking for convergence against specified criteria. Due to fundamental differences in their representation and computation, distinct paths are often taken within the SCF cycle for Plane-Wave (PW) and Atomic-Centered (AC) basis sets.

## Key Components

- **Module:** `esl_scf_m`
    - **Description:** Orchestrates the Self-Consistent Field (SCF) iterative process to determine the ground state electronic density and energy.
- **Type:** `scf_t`
    - **Description:** A derived type that encapsulates all data and configurations required for an SCF calculation.
    - **Components:**
        - `tol_reldens (real(dp))`: The tolerance threshold for the relative difference between input and output densities, used as the convergence criterion.
        - `max_iter (integer(ip))`: The maximum number of SCF iterations allowed before the loop is terminated if convergence is not met.
        - `mixer (type(mixing_t))`: An object of `mixing_t` (from `esl_mixing_m`) used to mix densities from successive iterations, which helps in stabilizing and accelerating convergence.
        - `rho_in (type(density_t))`: The input electronic density for a given SCF iteration.
        - `rho_out (type(density_t))`: The output electronic density calculated within an SCF iteration based on `rho_in`.
        - `H (type(hamiltonian_t))`: The Hamiltonian object, which includes the potential and is used to solve for electronic states.
    - **Procedures (public):**
        - `init`: Initializes the `scf_t` object and its components.
        - `mix`: Performs the density mixing step.
        - `loop`: Executes the main SCF iterative cycle.

### Subroutine/Function: `init(this, system, states)`

- **Description:** Sets up the `scf_t` object for an SCF calculation. Key actions include:
    1.  Reading SCF parameters from FDF input: `SCFTolerance` (for `this%tol_reldens`) and `SCFMaxIterations` (for `this%max_iter`).
    2.  Initializing the density mixer: `call this%mixer%init()`.
    3.  Performing basis-specific SCF initializations (currently placeholder comments exist for PW and AC specific initializations within this routine, but core AC setup like sparse patterns and overlap seem to be handled externally or assumed prior).
    4.  Initializing the Hamiltonian: `call this%H%init(system%basis, system%geo, states, periodic=.true.)`. This also initializes the `potential_t` member of `this%H`. For AC basis, this step typically includes computing fixed, non-SCF parts of the Hamiltonian.
    5.  Initializing the input (`this%rho_in`) and output (`this%rho_out`) density objects (`density_t`). This is done after Hamiltonian initialization as AC density matrices might depend on data (like sparse patterns) prepared during Hamiltonian setup.
    6.  If `WITH_FLOOK` is enabled, registers key SCF parameters (mixer alpha, max iterations, tolerance) to a Lua dictionary for scripting access.
- **Arguments:**
    - `this (class(scf_t))`: The SCF object to be initialized.
    - `system (type(system_t), intent(inout))`: The overall system configuration.
    - `states (type(states_t), intent(in))`: Initial electronic states information.
- **Returns:** Not applicable (Subroutine).

### Subroutine/Function: `loop(this, elsi, system, states, smear)`

- **Description:** This is the core routine that executes the iterative SCF cycle.
    1.  An initial guess for the input density `this%rho_in` is made using `this%rho_in%guess(system)`.
    2.  The potential part of `this%H` is calculated based on this initial guess. For AC basis, `this%H%ac%calculate_H0` is called (note: also called in `H%init`). For PW basis, `states` are randomized.
    3.  The main SCF loop (`loop_scf`) begins, iterating from `iter = 1` to `this%max_iter`:
        *   **Basis-Specific SCF Step:** Calls either `SCF_AC()` or `SCF_PW()` (internal subroutines) to perform the bulk of the work for one SCF iteration.
        *   **Energy Calculation:** `call system%energy%calculate()` computes various energy components based on the current electronic structure.
        *   **Convergence Check:** The residue `res = this%rho_in%residue(system%basis, this%rho_out, states)` is calculated, quantifying the difference between the input density of the current iteration and the output density just computed.
        *   **Logging:** Iteration number, total energy, Fermi level, and residue are logged using YAML.
        *   **Lua Hook:** If `WITH_FLOOK` is enabled, `call flook_if_call(LUA, LUA_SCF_LOOP)` allows Lua scripts to interact mid-loop.
        *   **Exit Condition:** If `res <= this%tol_reldens`, the SCF cycle is considered converged, and the loop is exited.
        *   **Density Mixing:** If not converged, `call this%mix(system%basis, this%rho_in, this%rho_out)` is performed. This updates `this%rho_in` with a mixed version of the old `rho_in` and the new `rho_out`, preparing it for the next iteration.
    4.  After the loop terminates (either by convergence or max iterations), final system energies are displayed using `system%energy%display()`.
- **Internal Subroutines:**
    - `SCF_AC()`: Details one iteration for an Atomic-Centered basis. Involves: calculating fixed energy terms, setting up `H0` in `this%H%ac`, calculating density on grid from `rho_in%ac%DM`, computing potentials (Hartree, XC, external) and adding them to `this%H%ac`, solving eigensystem via `this%H%eigensolver`, updating occupations and Fermi level via `smear%calc_fermi_occ`, and finally calculating the new `this%rho_out%ac%DM` (density matrix) from the updated states.
    - `SCF_PW()`: Details one iteration for a Plane-Wave basis. Involves: solving eigensystem via `this%H%eigensolver`, updating occupations and Fermi level via `smear%calc_fermi_occ`, calculating the new `this%rho_out%pw%rho` (density on grid) from updated states, and then recalculating potentials based on this new density.
- **Arguments (for `loop`):**
    - `this (class(scf_t), intent(inout))`: The SCF control object.
    - `elsi (type(elsi_t), intent(inout))`: ELSI library interface object, used by eigensolvers and smearing.
    - `system (type(system_t), intent(inout))`: The system object, its energy components are updated.
    - `states (type(states_t), intent(inout))`: Electronic states, updated by the eigensolver.
    - `smear (type(smear_t), intent(inout))`: Smearing object, used to calculate Fermi level and occupations.
- **Returns:** Not applicable (Subroutine).

### Subroutine/Function: `mix(this, basis, in, out)`

- **Description:** Applies a mixing scheme to the input (`in`) and output (`out`) densities of an SCF iteration to produce the input density for the next iteration (which is stored back into `in`). This is crucial for achieving stable convergence. It uses the `linear` mixing method from `this%mixer`. The actual data mixed depends on the basis type: for PW, it's the real-space charge density `rho`; for AC, it's the density matrix elements `DM%M`.
- **Arguments:**
    - `this (class(scf_t))`: The SCF object, to access `this%mixer`.
    - `basis (type(basis_t), intent(in))`: Basis set information, to determine data layout (PW grid size or AC DM non-zero elements).
    - `in (type(density_t), intent(inout))`: On input, the density from the previous iteration (or initial guess). On output, the newly mixed density to be used for the next iteration.
    - `out (type(density_t), intent(inout))`: The density computed in the current iteration before mixing.
- **Returns:** Not applicable (Subroutine).

## Important Variables/Constants

- **Within `scf_t` type:**
    - `tol_reldens (real(dp))`: Defines the convergence threshold for the SCF cycle based on density change.
    - `max_iter (integer(ip))`: Sets the maximum number of iterations for the SCF loop.

## Usage Examples

```fortran
! TODO: Add usage example
! Conceptual usage within a larger program (e.g., esl_main):
!
! use esl_scf_m
! ! ... other necessary modules ...
! type(scf_t) :: scf_calc
! type(system_t) :: current_system
! type(states_t) :: current_states
! type(elsi_t) :: elsi_interface
! type(smear_t) :: smearing_params
!
! ! ... (initialize system, states, elsi, smear) ...
!
! call scf_calc%init(current_system, current_states)
! call scf_calc%loop(elsi_interface, current_system, current_states, smearing_params)
!
! ! SCF calculation is now complete (or max_iter reached).
! ! Results are in current_system%energy, current_states, etc.
```

## Dependencies and Interactions

- **Core Modules Interacted With:**
    - `esl_density_m`: `rho_in` and `rho_out` are of `density_t`. Methods like `guess`, `calculate`, `calculate_density_matrix`, and `residue` are heavily used.
    - `esl_hamiltonian_m`: `this%H` is of `hamiltonian_t`. Its `init` and `eigensolver` methods are central. Components like `H%potential` and basis-specific parts (`H%ac`, `H%pw`) are directly accessed.
    - `esl_potential_m`: Used via `this%H%potential` to calculate various potential terms.
    - `esl_states_m`: `states_t` objects are passed in and modified by eigensolvers.
    - `esl_system_m`: `system_t` provides context (basis, geometry) and stores results like energy.
    - `esl_mixing_m`: `this%mixer` (of `mixing_t`) is used for density mixing.
    - `esl_smear_m`: `smear_t` is used for Fermi-level determination and occupation calculation.
    - `esl_basis_m`: Provides `basis_t` and type constants (`PLANEWAVES`, `ATOMCENTERED`) for dispatching logic.
- **Input/Output & Control:**
    - `fdf`: Reads SCF control parameters (`SCFTolerance`, `SCFMaxIterations`).
    - `yaml_output`: Used for structured logging of SCF progress (iteration details, energies, residue).
- **Optional Features:**
    - `esl_flook_...` modules: If `WITH_FLOOK` is defined, enables Lua scripting hooks within the SCF loop.
- **Atomic-Centered Specific Modules:**
    - `esl_sparse_pattern_m`, `esl_create_sparse_pattern_ac_m`, `esl_overlap_matrix_ac_m`, `esl_density_matrix_ac_m`: These modules are relevant for the setup and manipulation of quantities specific to atomic-centered basis sets, particularly the density matrix and Hamiltonian construction.
- **Overall Flow:**
    The `esl_scf_m` module is typically driven by a higher-level routine (e.g., `esl_main`). After `init` is called, the `loop` method takes control. Inside the loop, it iteratively updates the potential (based on `rho_in`), solves the eigenvalue problem (updating `states` and `rho_out`), checks for convergence, and mixes densities until `tol_reldens` is met or `max_iter` is exceeded.
