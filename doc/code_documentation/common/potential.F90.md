## Overview

The `esl_potential_m` module is responsible for managing and computing the various components of the effective potential experienced by electrons in an electronic structure calculation, particularly within the framework of Density Functional Theory (DFT). It defines a derived type, `potential_t`, which stores these potential components—Hartree potential, external local potential (typically from atomic nuclei/pseudopotentials), and the exchange-correlation (XC) potential—on a real-space grid. The module includes functionalities to initialize these potentials, compute them based on the electron density, and manage sub-components like a Poisson solver and an XC functional handler.

## Key Components

- **Module:** `esl_potential_m`
    - **Description:** Defines and manages different components of the electronic potential (Hartree, external, XC) on a real-space grid.
- **Type:** `potential_t`
    - **Description:** A derived type that encapsulates various potential arrays and related calculation utilities.
    - **Components:**
        - `np (integer)`: The number of points in the real-space grid on which the potentials are defined.
        - `hartree (real(dp), allocatable :: hartree(:))`: A 1D array storing the Hartree potential at each grid point. This potential arises from the classical electrostatic repulsion between electrons.
        - `external (real(dp), allocatable :: external(:))`: A 1D array storing the total external local potential at each grid point. This typically includes the local part of ionic pseudopotentials from all atoms.
        - `vxc (real(dp), allocatable :: vxc(:))`: A 1D array storing the exchange-correlation potential at each grid point, derived from the chosen XC functional.
        - `ionicOffset (real(dp))`: An offset value related to the external potential, possibly used in the Poisson solver for numerical stability or handling specific boundary conditions (e.g., Ewald sum related offsets).
        - `psolver (type(psolver_t))`: An object of type `psolver_t` (from `esl_psolver_m`), used to solve the Poisson equation to obtain the Hartree potential from the electron density.
        - `xc (type(xc_t))`: An object of type `xc_t` (from `esl_xc_m`), used to calculate the exchange-correlation potential and energy from the electron density.
    - **Procedures (public, type-bound):**
        - `init`: Initializes the potential object, its arrays, and its `psolver` and `xc` components. Also pre-computes the `external` local potential.
        - `calculate`: Computes the `hartree` and `vxc` potentials based on a given electron `density`.
        - `compute_ext_loc`: Specifically computes the `external` local potential by summing contributions from all atoms.
    - **Finalizer (type-bound):**
        - `finalizer`: Deallocates the `hartree`, `external`, and `vxc` arrays.

### Type-Bound Procedures

#### `init(pot, basis, states, geo, periodic)`
- **Description:** Initializes the `potential_t` object `pot`.
    1.  Sets `pot%np` (number of grid points) from `basis%grid%np`.
    2.  Allocates `pot%hartree`, `pot%external`, and `pot%vxc` arrays to size `pot%np` and initializes them to zero.
    3.  Determines a `geocode` ('P' for periodic, 'F' for free/isolated based on the `periodic` flag) for configuring the Poisson solver.
    4.  Sets `pot%ionicOffset = 0.0_dp`.
    5.  Initializes the embedded Poisson solver: `call pot%psolver%init(...)`, providing grid dimensions and spacings.
    6.  Initializes the embedded XC handler: `call pot%xc%init(geo, basis%grid)`.
    7.  Calls `pot%compute_ext_loc(basis, geo)` to pre-calculate and store the total external local potential from all atoms.
- **Arguments:**
    - `pot (class(potential_t), intent(inout))`: The potential object to be initialized.
    - `basis (class(basis_base_t), intent(in))`: Provides access to grid information via `basis%grid`.
    - `states (type(states_t), intent(in))`: Electronic states information (passed for interface consistency, not directly used in this specific `init` body).
    - `geo (type(geometry_t), intent(in))`: System geometry, used for `xc%init` and `compute_ext_loc`.
    - `periodic (logical, intent(in))`: Flag indicating if the system has periodic boundary conditions, used for `psolver` setup.

#### `calculate(pot, density, energy)`
- **Description:** Computes the Hartree potential and the exchange-correlation potential based on the provided electron `density`.
    1.  The input `density` array is first copied into `pot%hartree`. This array then serves as the source term (charge density) for the Poisson solver.
    2.  `call pot%psolver%hartree_potential(pot%hartree, pot%external, pot%ionicOffset, energy%hartree)`: This invokes the Poisson solver. It solves for the Hartree potential (which overwrites the density in `pot%hartree`) and calculates the Hartree energy, adding it to `energy%hartree`. The `pot%external` and `pot%ionicOffset` arguments might be used by the solver for specific configurations or corrections.
    3.  `call pot%xc%calculate(density, energy, pot%vxc)`: This calls the XC handler to compute the exchange-correlation potential (stored in `pot%vxc`) and the XC energy components (which are added to `energy%exchange` and `energy%correlation` within the `xc%calculate` method).
- **Arguments:**
    - `pot (class(potential_t), intent(inout))`: The potential object. Its `hartree` and `vxc` arrays are updated with the calculated potentials.
    - `density (real(dp), intent(in) :: density(:))`: The electron density on the real-space grid.
    - `energy (type(energy_t), intent(inout))`: The energy object where the calculated Hartree and XC energy contributions are accumulated.

#### `compute_ext_loc(pot, basis, geo)`
- **Description:** Calculates the total external local potential, `pot%external`, on the grid. This potential is typically the sum of the local parts of the ionic pseudopotentials from all atoms in the system.
    1.  Initializes `pot%external` to zero.
    2.  It iterates through each atom (`iat`) defined in `geo`.
    3.  For each atom, it retrieves its species information (`geo%species(is)%vlocal`, which is a tabulated radial local pseudopotential).
    4.  It calls `basis%grid%radial_function` to evaluate this radial pseudopotential on the grid, centered at the atom's position (`geo%xyz(:,iat)`). The result is stored in a temporary array `extloc`.
    5.  The values in `extloc` are then added to the corresponding points in `pot%external`.
- **Arguments:**
    - `pot (class(potential_t), intent(inout))`: The potential object; its `external` array is computed and stored.
    - `basis (class(basis_base_t), intent(in))`: Provides the grid definition via `basis%grid` and its `radial_function` method.
    - `geo (type(geometry_t), intent(in))`: Provides atomic species data (including `vlocal`) and atomic positions.

#### `finalizer(pot)`
- **Description:** This is a type-bound `final` subroutine, automatically called when a `potential_t` object is deallocated. It ensures that the dynamically allocated arrays `pot%hartree`, `pot%external`, and `pot%vxc` are deallocated if they were previously allocated, preventing memory leaks.
- **Arguments:**
    - `pot (type(potential_t))`: The `potential_t` object being finalized.

## Important Variables/Constants
- **Arrays:** `hartree(:)`, `external(:)`, `vxc(:)` are the core data stored, representing the different potential components on the grid.
- **Objects:** `psolver (type(psolver_t))` and `xc (type(xc_t))` are important embedded objects that provide the machinery for specific calculations.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage within an SCF iteration:
!
! use esl_potential_m
! use esl_density_m
! use esl_energy_m
! ! ... other necessary modules ...
!
! type(potential_t) :: sys_potential
! type(density_t) :: current_density  ! Assumed to have current_density%rho populated
! type(energy_t) :: current_energies
! type(basis_base_t) :: sys_basis_base ! Initialized with grid
! type(states_t) :: sys_states
! type(geometry_t) :: sys_geometry
! logical :: is_periodic
!
! ! Initialization (typically once before SCF loop)
! ! ... (initialize sys_basis_base, sys_states, sys_geometry, is_periodic) ...
! call sys_potential%init(sys_basis_base, sys_states, sys_geometry, is_periodic)
!
! ! Inside SCF loop:
! ! ... (current_density%rho is updated) ...
! ! ... (current_energies is initialized or carries over relevant parts) ...
! call sys_potential%calculate(current_density%rho, current_energies)
!
! ! Now sys_potential%hartree, sys_potential%vxc, and sys_potential%external
! ! can be used to construct the Kohn-Sham Hamiltonian.
! ! current_energies%hartree, current_energies%exchange, current_energies%correlation are updated.
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec`: For precision kinds `dp` (double precision) and `ip` (integer).
    - `esl_basis_base_m`: Provides `basis_base_t`, which is used to access grid information (`basis%grid`).
    - `esl_energy_m`: Provides `energy_t` for accumulating energy contributions (Hartree, XC energies).
    - `esl_geometry_m`: Provides `geometry_t` for atomic species information (local pseudopotentials `vlocal`) and positions.
    - `esl_grid_m`: The `grid_t` type (accessed via `basis%grid`) and its `radial_function` method are used in `compute_ext_loc`.
    - `esl_psolver_m`: Provides `psolver_t`, an instance of which (`pot%psolver`) is used to calculate the Hartree potential.
    - `esl_states_m`: Provides `states_t`, used in the interface of `init` (though not directly in its body).
    - `esl_xc_m`: Provides `xc_t`, an instance of which (`pot%xc`) is used to calculate the exchange-correlation potential and energy.
    - `esl_constants_m`: Potentially used indirectly via calls to methods of other ESL modules (e.g., `radial_function` might use `PI`).
- **External Libraries:**
    - **FFTW:** The Fast Fourier Transform library is likely an indirect dependency via `esl_psolver_m`, as Poisson solvers on a grid typically use FFTs.
- **Interactions with other components:**
    - **`esl_hamiltonian_m`:** The `potential_t` object is a critical input for constructing the Kohn-Sham Hamiltonian. The `hamiltonian_t` type (from `esl_hamiltonian_m`) usually holds a `potential_t` object. The sum of `pot%hartree`, `pot%external`, and `pot%vxc` forms the effective Kohn-Sham potential.
    - **`esl_scf_m` (SCF Cycle):** The `calculate` method of `potential_t` is called in each iteration of the Self-Consistent Field (SCF) loop. It takes the current electron `density` as input and computes the new potentials and their corresponding energy contributions.
    - **`esl_grid_m`:** All potential components are defined on a `grid_t` (from `esl_grid_m`). The `compute_ext_loc` method directly uses `grid_t%radial_function` to map atomic local pseudopotentials onto this grid.
    - **`esl_species_m`:** Information from `species_t` objects (specifically `species%vlocal`, the radial local pseudopotential) is used in `compute_ext_loc`.
    - **Memory Management:** The `finalizer` ensures proper deallocation of the large potential arrays.
