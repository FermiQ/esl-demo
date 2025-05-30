## Overview

The `esl_system_m` module defines and manages the `system_t` derived type, which serves as the primary top-level container for all data pertinent to an electronic structure calculation within this software package. As described in its comments, `system_t` "encompass[es] everything related to a single DFT calculation" and is considered the "highest level type which is necessary for performing any DFT calculation." It aggregates fundamental aspects of the simulated system, including its basis set representation, various energy components, atomic geometry (positions and cell parameters), forces on atoms, and ion-ion interaction energies.

## Key Components

- **Module:** `esl_system_m`
    - **Description:** Provides the `system_t` derived type, which acts as a central repository for all information defining a computational electronic structure problem.
- **Type:** `system_t`
    - **Description:** A comprehensive derived type that groups together various constituent objects representing different facets of the system being simulated.
    - **Components:**
        - `basis (type(basis_t))`: An object of `basis_t` (from `esl_basis_m`) that stores all information related to the chosen basis set (e.g., Plane-Wave or Atomic-Centered).
        - `energy (type(energy_t))`: An object of `energy_t` (from `esl_energy_m`) that holds the values of various energy components of the system (e.g., total energy, Hartree energy, exchange-correlation energy, kinetic energy, Fermi energy, ion-ion energy).
        - `geo (type(geometry_t))`: An object of `geometry_t` (from `esl_geometry_m`) that describes the system's geometry, including atomic species, coordinates, and unit cell vectors if applicable.
        - `force (type(force_t))`: An object of `force_t` (from `esl_force_m`) used to store the forces acting on each atom in the system.
        - `ion_inter (type(ion_interaction_t))`: An object of `ion_interaction_t` (from `esl_ion_interaction_m`) responsible for calculating the interaction energy between atomic nuclei (ions) and the corresponding forces.
        - `nElectrons (real(dp))`: A scalar variable intended to store the total number of electrons in the system. It is initialized to `0.0_dp`, suggesting it's populated during a more detailed setup phase (e.g., from species information or explicitly set).
    - **Procedures (public):**
        - `init`: Initializes the `system_t` object and all its member components.
        - `update`: Updates certain system properties, particularly ion-ion interactions, typically after a change in atomic geometry (e.g., in a molecular dynamics step).
        - `summary`: Prints a concise summary of the system's configuration (geometry and basis) in YAML format.

### Subroutine/Function: `init(this)`

- **Description:** This subroutine initializes a `system_t` object by sequentially calling the `init` methods of its constituent member objects: `this%energy`, `this%geo`, `this%basis` (which requires `this%geo` as input), `this%force` (which requires the number of atoms from `this%geo%n_atoms`), and `this%ion_inter`. The call to `this%energy%init()` appears twice; this might be an oversight or serve a specific reset purpose. After individual initializations, `this%summary()` is called to output a summary of the newly configured system.
If the code is compiled `WITH_FLOOK` (Lua scripting enabled), this routine also registers a wide array of system variables—including atomic coordinates, cell vectors, various force components, Ewald parameter for ion-ion interaction, and multiple energy terms—into a Lua dictionary. This makes these critical system parameters accessible and potentially modifiable from Lua scripts.
- **Arguments:**
    - `this (class(system_t), intent(inout))`: The `system_t` object to be initialized.
- **Returns:** Not applicable (Subroutine).

### Subroutine/Function: `update(this, periodic)`

- **Description:** This subroutine is designed to be called when the system's ionic configuration changes (e.g., after an MD step) but before a new electronic SCF calculation is performed. It recalculates quantities that depend directly on the ionic geometry but not on the electron density from the SCF cycle. Specifically, it updates the ion-ion interaction energy (`this%energy%ionion`) and the force component arising from ion-ion interactions (`this%force%ionion`). It does this by calling either `this%ion_inter%calculate_periodic` or `this%ion_inter%calculate_isolated`, depending on the boolean `periodic` flag.
- **Arguments:**
    - `this (class(system_t), intent(inout))`: The `system_t` object to be updated.
    - `periodic (logical, intent(in))`: A flag that determines whether to use calculation methods for periodic systems (e.g., Ewald summation) or isolated systems.
- **Returns:** Not applicable (Subroutine).

### Subroutine/Function: `summary(sys)`

- **Description:** Outputs a summary of the system's setup to the standard output stream using YAML format. It creates a top-level YAML mapping called "System" and then calls the `summary` methods of the `sys%geo` (geometry) and `sys%basis` (basis set) components to include their respective details.
- **Arguments:**
    - `sys (class(system_t))`: The `system_t` object for which the summary is to be generated.
- **Returns:** Not applicable (Subroutine).

## Important Variables/Constants

- **`nElectrons (real(dp))`**: A member of the `system_t` type, designated for storing the total number of electrons in the simulated system. Its actual value would be determined during the setup based on atomic species or explicit input.

## Usage Examples

```fortran
! TODO: Add usage example
! Conceptual usage:
!
! use esl_system_m
! type(system_t) :: my_computational_system
! logical :: is_periodic_system
!
! ! Initialize the entire system
! call my_computational_system%init()
!
! ! ... (perform some operations, e.g., an SCF calculation) ...
!
! ! If atomic positions change (e.g., MD step)
! is_periodic_system = .true. ! or .false.
! ! ... (update my_computational_system%geo%xyz) ...
! call my_computational_system%update(periodic=is_periodic_system)
!
! ! ... (proceed with next SCF or other calculations) ...
!
! ! Print a summary at any point
! call my_computational_system%summary()
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec`: Provides precision kinds (`dp`, `ip`).
    - `esl_basis_m`: Provides `basis_t` and its methods.
    - `esl_energy_m`: Provides `energy_t` and its methods.
    - `esl_force_m`: Provides `force_t` and its methods.
    - `esl_geometry_m`: Provides `geometry_t` and its methods.
    - `esl_ion_interaction_m`: Provides `ion_interaction_t` and its calculation routines.
    - `esl_species_m`: Used by `esl_geometry_m` for defining atomic species.
    - `yaml_output`: Used by the `summary` method for structured output.
    - `dictionary`, `esl_dict_m` (Conditional, if `WITH_FLOOK` is defined): Used for integrating with Lua scripting by exposing system variables.
- **External Libraries:**
    - `Lua` (Conditional, if `WITH_FLOOK` is defined): If enabled, the Lua library would be linked for scripting capabilities.
- **Interactions with other components:**
    - **Central Data Hub:** `system_t` acts as the central data structure passed around to various computational modules (e.g., `esl_scf_m`, `esl_density_m`, `esl_hamiltonian_m`, `esl_main`). These modules read from and write to the components of `system_t` (like `system%energy` or `system%geo`).
    - **Initialization Driver:** The `init` method of `system_t` drives the initialization of all its core components.
    - **Inter-Step Updates:** The `update` method is crucial for simulations involving multiple steps (like Molecular Dynamics or geometry optimizations), where it's called by modules like `esl_next_step_m` to refresh geometry-dependent but non-SCF quantities.
    - **Scripting Interface:** When Lua is enabled, `system_t%init` makes many fundamental physical parameters of the system available to Lua scripts, allowing for advanced control, analysis, or modification of the simulation's behavior.
    - The comment "TODO decide whether the system should be inherited for the LO/PW case" hints at potential future design considerations for handling systems with mixed or more complex basis set strategies (e.g., combined Localized Orbitals and Plane Waves).
