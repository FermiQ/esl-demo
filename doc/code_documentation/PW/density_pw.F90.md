## Overview

The `esl_density_pw_m` module defines the `density_pw_t` derived type, which is responsible for managing the electron density within a plane-wave (PW) computational scheme. This type extends `density_base_t`, inheriting its primary component `rho`, an allocatable array that stores the electron density on a real-space grid.

The `density_pw_t` type specifically includes a pointer to a `basis_pw_t` object, which provides the necessary context for plane-wave operations, such as G-vector mappings and FFT grid information. The module provides methods for:
-   Initializing the density object and linking it to a specific PW basis.
-   Generating an initial guess for the real-space density `rho` by summing tabulated atomic valence charge densities.
-   Calculating the real-space density `rho` from the electronic states (wavefunctions represented by their plane-wave coefficients, along with occupation numbers and k-point weights). This involves transforming wavefunction coefficients from reciprocal space to the real-space grid.
-   Computing the residual (difference) between two real-space density distributions.
-   Utility routines to get or set the real-space density `rho` from/to an external array.

## Key Components

- **Module:** `esl_density_pw_m`
    - **Description:** Manages the electron density for plane-wave (PW) basis calculations, including its representation on a real-space grid and transformations from PW coefficients.

- **Type:** `density_pw_t`, extends `density_base_t`
    - **Description:** A derived type for handling electron density in a PW framework.
    - **Inherited Components (from `density_base_t`):**
        - `rho (real(dp), allocatable :: rho(:))`: Stores the electron density values on a real-space grid.
    - **Specific Components:**
        - `np (integer)`: The number of points in the real-space grid (copied from `basis%grid%np` for convenience).
        - `pw (type(basis_pw_t), pointer :: pw => null())`: A pointer to the associated `basis_pw_t` object, which defines the plane-wave basis set and its properties (e.g., G-vectors, FFT grid).
    - **Procedures (Public, Type-Bound):**
        - `init(this, basis)`: Initializes the density object and associates it with a PW basis.
        - `guess(this, geo)`: Generates an initial guess for `this%rho` by summing atomic densities.
        - `calculate(this, states)`: Computes `this%rho` from electronic states (PW coefficients, occupations).
        - `residue(this, grid, nel, other)`: Calculates the integrated absolute difference between `this%rho` and `other%rho`.
        - `get_den(this, rho_out)`: Copies `this%rho` to an external array `rho_out`.
        - `set_den(this, rho_in)`: Sets `this%rho` from an external array `rho_in`.
    - **Finalizer (Type-Bound):**
        - `finalizer(this)`: Deallocates `this%rho` and nullifies the `this%pw` pointer.

### Type-Bound Procedures

#### `init(this, basis)`
- **Description:** Initializes a `density_pw_t` object. It allocates the real-space density array `this%rho` to match the size of the grid provided by `basis%grid%np`. It also stores this size in `this%np` and associates the pointer `this%pw` with the input `basis` object.
- **Arguments:**
    - `this (class(density_pw_t), intent(inout))`: The density object to initialize.
    - `basis (type(basis_pw_t), target, intent(in))`: The plane-wave basis set associated with this density.

#### `guess(this, geo)`
- **Description:** Generates an initial guess for the real-space electron density `this%rho`.
    1.  `this%rho` is initialized to all zeros.
    2.  It then iterates through each atom `iat` in the system geometry `geo`.
    3.  For each atom, it retrieves its tabulated atomic valence charge density (a `pspiof_meshfunc_t` from `geo%species(is)%rho`).
    4.  `this%pw%grid%radial_function` is called to evaluate this radial atomic density on the real-space grid, centered at the atom's position, storing the result in a temporary array `atomicden`.
    5.  The `atomicden` is then added to `this%rho`.
    The norm of each individual atomic density contribution is printed to YAML output during this process.
- **Arguments:**
    - `this (class(density_pw_t), intent(inout))`: The density object; `this%rho` is populated.
    - `geo (type(geometry_t), intent(in))`: The system geometry, providing atomic species (and their reference charge densities) and positions.

#### `calculate(this, states)`
- **Description:** Computes the total real-space electron density `this%rho` from the electronic states (`states_t` object).
    1.  `this%rho` is initialized to all zeros.
    2.  It iterates over all k-points (`ik`), spin components (`isp`), and electronic states/bands (`ist`) defined in the `states` object.
    3.  For each state $(ist, isp, ik)$:
        *   It calls `pw2grid` (from `esl_utils_pw_m`) to transform the plane-wave coefficients of the wavefunction (`states%states(ist,isp,ik)%zcoef`) from reciprocal space to its real-space representation (`coef_rs`) on the grid defined by `this%pw%grid`. This transformation uses G-vector information from `this%pw%gmap` and `this%pw%ndims`.
        *   The contribution of this state to the density at each grid point `ip` is calculated as `weight * abs(coef_rs(ip))**2`, where `weight = states%occ_numbers(ist,isp,ik) * states%k_weights(ik)`.
        *   This contribution is accumulated into `this%rho(ip)`.
    A `TODO` in the source code notes that `this%pw%gmap` should ideally be k-point specific if more than the Gamma point is used.
- **Arguments:**
    - `this (class(density_pw_t), intent(inout))`: The density object; `this%rho` is computed.
    - `states (type(states_t), intent(in))`: The electronic states, providing PW coefficients, occupation numbers, eigenvalues, and k-point information.

#### `residue(this, grid, nel, other) result(res)`
- **Description:** Calculates a measure of the difference between two real-space density distributions, `this%rho` and `other%rho`. It computes the integral of the absolute difference: `res = integral |other%rho(r) - this%rho(r)| dr`. The integration is performed numerically using `grid%integrate`. If `nel` (total number of electrons) is greater than zero, the result is normalized by `nel`.
- **Arguments:**
    - `this (class(density_pw_t), intent(in))`: The first density object.
    - `grid (type(grid_t), intent(in))`: The grid object used for integration.
    - `nel (real(dp), intent(in))`: The total number of electrons, for normalization.
    - `other (type(density_pw_t), intent(in))`: The second density object for comparison.
- **Returns:** `res (real(dp))`: The calculated (and optionally normalized) residue.

#### `get_den(this, rho_out)`
- **Description:** Copies the real-space density `this%rho` into an external array `rho_out`.
- **Arguments:**
    - `this (class(density_pw_t), intent(in))`.
    - `rho_out (real(dp), intent(out) :: rho_out(:))`: The output array to receive the density. Must be allocated to the correct size (`this%np`).

#### `set_den(this, rho_in)`
- **Description:** Sets the real-space density `this%rho` from an external array `rho_in`.
- **Arguments:**
    - `this (class(density_pw_t), intent(inout))`.
    - `rho_in (real(dp), intent(in) :: rho_in(:))`: The input array providing the density. Must have size `this%np`.

#### `finalizer(this)`
- **Description:** Type-bound final procedure. Deallocates `this%rho` if it is allocated and nullifies the `this%pw` pointer.
- **Arguments:** `this (type(density_pw_t), intent(inout))`.

## Important Variables/Constants
- **`rho(:)`**: (Inherited) The array storing the real-space electron density.
- **`pw` (pointer to `basis_pw_t`)**: Links this density object to its corresponding plane-wave basis definition, crucial for PW-specific operations.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage:
!
! use esl_density_pw_m
! use esl_basis_pw_m
! use esl_geometry_m
! use esl_states_m
!
! type(density_pw_t) :: pw_density
! type(basis_pw_t) :: pw_basis
! type(geometry_t) :: system_geometry
! type(states_t) :: current_electron_states
!
! ! ... (Initialize pw_basis, system_geometry, current_electron_states) ...
!
! ! Initialize density object
! call pw_density%init(pw_basis)
!
! ! Generate initial guess for density
! call pw_density%guess(system_geometry)
!
! ! In SCF loop, after states are updated:
! call pw_density%calculate(current_electron_states)
!
! ! ... (Use pw_density%rho for potential calculation, check convergence with residue) ...
!
! ! pw_density is automatically finalized when it goes out of scope.
```

## Dependencies and Interactions

- **Base Type:** `esl_density_base_m`: `density_pw_t` extends `density_base_t`.
- **Core ESL Modules:**
    - `prec`: For precision kinds (`dp`).
    - `esl_basis_pw_m`: Provides `basis_pw_t`, which is essential for PW-specific operations like the `pw2grid` transformation in `calculate` and for accessing grid information.
    - `esl_geometry_m`: Provides `geometry_t`, used in `guess` to access atomic species' reference charge densities and positions.
    - `esl_grid_m`: The `grid_t` type (accessed via `this%pw%grid`) and its methods (e.g., `radial_function`, `integrate`) are used.
    - `esl_states_m`: Provides `states_t`, which is the input for the `calculate` method (containing PW coefficients, occupations, etc.).
    - `esl_constants_m`: Provides `PI`, used in `guess` (though its direct use in the snippet is not apparent, it's often needed for normalization or spherical averaging in atomic density generation).
    - `yaml_output`: Used in `guess` for logging.
- **PW Specific Utility Modules:**
    - `esl_utils_pw_m`: Provides `pw2grid` for transforming PW coefficients to the real-space grid.
- **Interactions with other components:**
    - **`esl_density_m` (Main Density Wrapper):** `density_pw_t` is the specialized type for PW calculations, likely held within an instance of `esl_density_m%density_t%pw`.
    - **SCF Cycle (e.g., `esl_scf_m`):**
        - `guess` is called to provide an initial real-space density.
        - `calculate` is called in each SCF iteration after new electronic states (PW coefficients and occupations) are determined, to compute the updated real-space density.
        - `residue` can be used to compare densities from successive iterations for SCF convergence checking.
        - The computed `this%rho` is then used by `esl_potential_m` to calculate the Hartree and XC potentials for the next SCF step.
    - **Potential Calculation (`esl_potential_m`):** The real-space density `this%rho` is the primary input for calculating the Hartree and exchange-correlation potentials.
    - **Basis Set (`esl_basis_pw_m`):** The `density_pw_t` is tightly coupled with a `basis_pw_t` object, which provides the G-vector definitions, FFT grid, and other parameters necessary for PW operations.
</tbody>
</table>
