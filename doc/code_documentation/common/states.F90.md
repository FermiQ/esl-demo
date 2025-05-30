## Overview

The `esl_states_m` module is dedicated to managing the electronic states of a quantum mechanical system. It defines derived types to store all relevant information about these states, including the wavefunction coefficients, occupation numbers, eigenvalues (energies), and details about spin polarization and k-point sampling.

The module introduces two primary types:
-   `wfn_t`: A simple container for the coefficients of a single electronic wavefunction. These coefficients can be real or complex.
-   `states_t`: A comprehensive container that holds arrays of `wfn_t` objects for all states, along with corresponding arrays for their occupation numbers and eigenvalues, distributed over spin components and k-points. It also stores metadata like the total number of electrons and the number of coefficients per wavefunction.

Functionalities include initializing the state ensemble (determining the number of states based on electron count, allocating memory), randomizing wavefunction coefficients for initial guesses, providing a summary of the state configuration, and ensuring proper deallocation of memory.

## Key Components

- **Module:** `esl_states_m`
    - **Description:** Defines data structures and routines for managing electronic states (wavefunctions, occupations, eigenvalues).

- **Type:** `wfn_t`
    - **Description:** A simple derived type to hold the expansion coefficients of a single wavefunction in a chosen basis set.
    - **Components:**
        - `dcoef (real(dp), allocatable :: dcoef(:))`: A 1D allocatable array to store real-valued wavefunction coefficients.
        - `zcoef (complex(dp), allocatable :: zcoef(:))`: A 1D allocatable array to store complex-valued wavefunction coefficients. Only one of `dcoef` or `zcoef` would typically be allocated for a given wavefunction.

- **Type:** `states_t`
    - **Description:** The main derived type for managing the complete set of electronic states in the system.
    - **Components:**
        - `nstates (integer)`: The number of electronic states (bands) being considered.
        - `nspin (integer)`: The number of spin components (typically 1 for non-spin-polarized or spin-restricted calculations, 2 for spin-polarized or spin-unrestricted).
        - `nkpt (integer)`: The number of k-points used in sampling the Brillouin zone (1 for Gamma-point calculations).
        - `nel (real(dp))`: The total number of electrons in the system.
        - `complex_states (logical)`: A flag that is `.true.` if the wavefunction coefficients are complex and `.false.` if they are real.
        - `ncoef (integer)`: The number of coefficients used to represent each wavefunction (e.g., the size of the basis set).
        - `states (type(wfn_t), allocatable :: states(:,:,:))`: A 3D allocatable array of `wfn_t` objects, dimensioned (`nstates, nspin, nkpt`), storing the coefficients for each electronic state.
        - `occ_numbers (real(dp), allocatable :: occ_numbers(:,:,:))`: A 3D allocatable array of the same dimensions, storing the occupation number (typically between 0 and 2/`nspin`) for each state.
        - `k_weights (real(dp), allocatable :: k_weights(:))`: A 1D allocatable array of dimension (`nkpt`), storing the weight of each k-point.
        - `eigenvalues (real(dp), allocatable :: eigenvalues(:,:,:))`: A 3D allocatable array of the same dimensions, storing the eigenvalue (energy) of each state.
    - **Procedures (public, type-bound):**
        - `init`: Initializes the `states_t` object, allocates its arrays, and sets initial occupations.
        - `summary`: Prints a summary of the states configuration (number of states, spins, k-points).
        - `randomize`: Fills wavefunction coefficients with random numbers and normalizes them, useful for an initial guess.
    - **Finalizer (type-bound):**
        - `finalizer`: Handles deallocation of all allocatable arrays within the `states_t` object, including the coefficient arrays within each `wfn_t`.

### Type-Bound Procedures of `states_t`

#### `init(this, geo, nspin, nkpt, complex_wfn, num_coeffs)`
- **Description:** Initializes a `states_t` object `this`.
    1.  Sets basic parameters: `this%nspin`, `this%nkpt`, `this%complex_states = complex_wfn`, `this%ncoef = num_coeffs`.
    2.  Determines the total number of electrons `this%nel` by calling `geo%electronic_charge()`.
    3.  Calculates the minimum number of states `this%nstates` required to hold `this%nel` electrons (considering `this%nspin`). It then adds any `extra_states` specified in the FDF input file (key `ExtraStates`, default 0).
    4.  Allocates the main arrays: `this%states`, `this%occ_numbers`, and `this%eigenvalues` to dimensions (`this%nstates, this%nspin, this%nkpt`).
    5.  Sets initial `this%occ_numbers`: It fills states with the maximum possible occupation for the given spin polarization (e.g., 2.0 for `nspin=1`, 1.0 for `nspin=2`) until all electrons (`this%nel`) are accounted for. Remaining states are given zero occupation.
    6.  Initializes `this%eigenvalues` to `0.0_dp`.
    7.  Allocates `this%k_weights` and initializes them to give equal weight to each k-point, summing to 1.0.
    8.  For each `wfn_t` object within `this%states(:,:,:)`, it allocates either its `dcoef` array (if `complex_wfn` is false) or its `zcoef` array (if `complex_wfn` is true) to the size `num_coeffs`.
- **Arguments:**
    - `this (class(states_t), intent(inout))`: The `states_t` object to initialize.
    - `geo (type(geometry_t), intent(inout))`: Geometry object, used to get the total electronic charge.
    - `nspin (integer, intent(in))`: Number of spin components.
    - `nkpt (integer, intent(in))`: Number of k-points.
    - `complex_wfn (logical, intent(in))`: Flag indicating if wavefunctions are complex-valued.
    - `num_coeffs (integer, intent(in))`: Number of coefficients per wavefunction (size of basis).

#### `randomize(this)`
- **Description:** Populates the wavefunction coefficients with pseudo-random numbers and then normalizes each state. This is often used to generate an initial guess for the wavefunctions at the start of an SCF calculation.
    - If `this%complex_states` is true, the `zcoef` array of each state is filled with complex numbers where the real and imaginary parts are `random_number() - 0.5`.
    - If `this%complex_states` is false, the `dcoef` array of each state is filled with real numbers `random_number() - 0.5`.
    - After populating with random numbers, each state (i.e., each `dcoef` or `zcoef` array) is normalized such that the sum of the absolute squares of its coefficients is 1.
- **Arguments:**
    - `this (class(states_t), intent(inout))`: The `states_t` object whose wavefunction coefficients are to be randomized.

#### `summary(this)`
- **Description:** Prints a brief summary of the electronic states configuration to standard output using YAML formatting. The summary includes the number of states, number of spin components, and number of k-points.
- **Arguments:**
    - `this (class(states_t), intent(inout))`: The `states_t` object to summarize. *(Note: `intent(in)` would be sufficient as `this` is not modified)*.

#### `finalizer(this)`
- **Description:** This is a type-bound `final` subroutine, automatically invoked when a `states_t` object is deallocated. It meticulously deallocates all dynamically allocated memory:
    - It iterates through each `wfn_t` object in the `this%states` array and deallocates its `dcoef` or `zcoef` array (whichever is allocated).
    - Then, it deallocates the main arrays: `this%states`, `this%occ_numbers`, `this%eigenvalues`, and `this%k_weights`.
- **Arguments:**
    - `this (type(states_t), intent(inout))`: The `states_t` object being finalized.

## Important Variables/Constants
- All components of `states_t` are essential for a complete description of the electronic state ensemble, including `nstates`, `nspin`, `nkpt`, `nel`, `complex_states`, `ncoef`, and the arrays `states` (containing `wfn_t`), `occ_numbers`, `k_weights`, and `eigenvalues`.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage:
!
! use esl_states_m
! use esl_geometry_m
! use prec, only: dp, ip
!
! type(states_t) :: electron_states
! type(geometry_t) :: system_geometry
! integer :: num_spin = 1, num_kpts = 1, num_basis_coeffs = 100
! logical :: use_complex_wfn = .false.
!
! ! ... (Initialize system_geometry, e.g., system_geometry%init() ) ...
!
! ! Initialize electronic states
! call electron_states%init(system_geometry, num_spin, num_kpts, use_complex_wfn, num_basis_coeffs)
!
! ! Randomize initial wavefunctions (optional, for starting SCF)
! call electron_states%randomize()
!
! ! Print summary of states
! call electron_states%summary()
!
! ! ... (electron_states would then be used in SCF, eigensolvers, density calculations) ...
!
! ! electron_states is automatically finalized when it goes out of scope.
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only : dp, ip`: Imports double precision (`dp`) and integer precision (`ip`) kind specifiers.
    - `fdf`: Used by `init` to get the `ExtraStates` parameter from an FDF input file.
    - `yaml_output`: Used by the `summary` procedure for YAML-formatted output.
    - `esl_geometry_m`: Provides `geometry_t`, which is used in `init` to determine the total number of electrons (`geo%electronic_charge()`).
- **External Libraries:**
    - This module uses the Fortran intrinsic `random_number` for the `randomize` procedure. It does not have direct dependencies on other external libraries.
- **Interactions with other components:**
    - **SCF Engine (e.g., `esl_scf_m`):** The `states_t` object is a central data structure in any Self-Consistent Field (SCF) calculation. It holds the wavefunctions that are iteratively refined, the eigenvalues that are solved for, and the occupation numbers that are determined at each step.
    - **Hamiltonian Solvers (e.g., within `esl_hamiltonian_m`):** Eigensolver routines take the current Hamiltonian and update the `eigenvalues` and wavefunction `coefficients` within the `states_t` object.
    - **Density Calculation (e.g., `esl_density_m`):** The electron density is computed from the `states%states` (wavefunction coefficients) and `states%occ_numbers`.
    - **Smearing Module (e.g., `esl_smear_m`):** For metallic systems, a smearing module would use the `eigenvalues` to calculate the `fermi_level` and update `occ_numbers`.
    - **Basis Set:** The `ncoef` member directly relates to the size of the basis set used to expand the Kohn-Sham orbitals (wavefunctions).
    - **Memory Management:** The `finalizer` is crucial for preventing memory leaks by ensuring all dynamically allocated components of `states_t` are properly deallocated.
</tbody>
</table>
