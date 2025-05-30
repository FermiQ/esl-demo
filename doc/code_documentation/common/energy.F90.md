## Overview

The `esl_energy_m` module provides a specialized derived type, `energy_t`, designed to store and manage the various energy components that are computed during an electronic structure calculation, particularly within a Density Functional Theory (DFT) framework. The module defines the structure for these energy terms, all of which are assumed to be in Hartree atomic units. It also provides utility subroutines to initialize these energy values, to calculate the total energy based on its constituent parts, and to display a formatted breakdown of all stored energy terms.

## Key Components

- **Module:** `esl_energy_m`
    - **Description:** Defines the `energy_t` data type for holding various system energy components and provides procedures for their management (initialization, calculation of total energy, display).
- **Type:** `energy_t`
    - **Description:** A derived type that serves as a container for different energy terms calculated in an electronic structure simulation. All energy components are public members of type `real(dp)` and are initialized to `0.0_dp`.
    - **Components (all `real(dp)`):**
        - `total`: The total energy of the system. The module comments emphasize that this value is updated specifically upon calling the `calculate` method.
        - `eigenvalues`: The sum of Kohn-Sham eigenvalues (often referred to as the band energy).
        - `hartree`: The classical electrostatic energy due to the electron-electron repulsion (Hartree energy).
        - `KB`: Energy contribution from Kleinman-Bylander non-local pseudopotentials.
        - `ionion`: The electrostatic interaction energy between atomic nuclei (ion-ion repulsion).
        - `extern`: Energy contribution from any external potential applied to the system.
        - `exchange`: The exchange part of the exchange-correlation energy functional.
        - `correlation`: The correlation part of the exchange-correlation energy functional.
        - `int_nvxc`: A term whose exact nature is marked with "TODO ?" in the source. It's often a correction term in DFT to avoid double-counting of electron-XC field interaction when using the sum of eigenvalues (e.g., related to `Integral n(r)*V_xc(r) dr`).
        - `kinetic`: The kinetic energy of the non-interacting Kohn-Sham electrons (T_s).
        - `entropy`: The electronic entropy term. The comment notes it as "unit-less", which might mean it's -S/k_B, or if it's -TS, then T is implicitly included. Typically, -TS is the energy contribution.
        - `fermi`: The Fermi level (chemical potential) of the electronic system.
    - **Procedures (public):**
        - `init`: Initializes all energy components to zero.
        - `calculate`: Computes the `total` energy from the other stored energy terms using a specific DFT formula.
        - `display`: Prints a formatted list of all energy components using YAML output.

### Subroutine/Function: `init(this)`
- **Description:** This subroutine sets all numerical members (energy components, entropy, Fermi level) of the `energy_t` object `this` to `0.0_dp`. It serves as a reset or initialization mechanism.
- **Arguments:**
    - `this (class(energy_t))`: The `energy_t` object to be initialized.
- **Returns:** Not applicable (Subroutine).

### Subroutine/Function: `calculate(this)`
- **Description:** Calculates the `this%total` energy based on the other energy components stored within the `energy_t` object. The formula used is:
  `this%total = this%ionion + this%eigenvalues + this%extern - this%hartree + this%exchange + this%correlation + this%kinetic - this%int_nvxc`
  This formula is characteristic of Kohn-Sham DFT, where the sum of eigenvalues includes the Hartree energy and interaction with the XC potential, necessitating corrections (subtracting `hartree` and `int_nvxc`) to avoid double counting.
- **Arguments:**
    - `this (class(energy_t))`: The `energy_t` object. Its `total` member will be updated. Other members are assumed to have been previously calculated and stored.
- **Returns:** Not applicable (Subroutine).

### Subroutine/Function: `display(this)`
- **Description:** Outputs a detailed breakdown of the energy components stored in the `energy_t` object `this`. The output is directed to the standard output stream and formatted using YAML conventions for clarity. Each energy term is listed with its corresponding value under a main "Energy" YAML mapping.
- **Arguments:**
    - `this (class(energy_t))`: The `energy_t` object whose components are to be displayed.
- **Returns:** Not applicable (Subroutine).

## Important Variables/Constants
- All named components within the `energy_t` type are significant as they represent distinct physical or computational contributions to the system's energy or electronic structure (e.g., `total`, `eigenvalues`, `hartree`, `kinetic`, `fermi`).

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage:
!
! use esl_energy_m
! type(energy_t) :: system_energies
!
! ! Initialize energy components
! call system_energies%init()
!
! ! ... (During an SCF calculation, populate individual energy terms) ...
! system_energies%eigenvalues = ... ! from eigensolver
! system_energies%hartree = ...     ! from Poisson solver
! system_energies%exchange = ...    ! from XC functional calculation
! system_energies%correlation = ... ! from XC functional calculation
! system_energies%kinetic = ...     ! from KS orbitals
! system_energies%ionion = ...      ! from ion positions
! system_energies%extern = ...      ! if external field present
! system_energies%int_nvxc = ...    ! E_xc - V_xc term or similar
! system_energies%fermi = ...       ! from occupation/smearing calculation
!
! ! Calculate the total energy
! call system_energies%calculate()
!
! ! Display the energy breakdown
! call system_energies%display()
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only : dp`: Imports the `dp` kind specifier for double precision real numbers, used for all energy terms.
    - `yaml_output`: This module is used by the `display` subroutine to produce YAML-formatted output of the energy components.
- **External Libraries:**
    - This module has no direct dependencies on external libraries.
- **Interactions with other components:**
    - **`esl_system_m`:** The `system_t` type (defined in `esl_system_m`) contains an instance of `energy_t` as one of its core members (usually accessed as `system%energy`). This makes `esl_energy_m` integral to the overall system description.
    - **SCF Cycle (e.g., `esl_scf_m`):** Throughout a Self-Consistent Field (SCF) calculation, various computational steps update the individual fields of the `energy_t` object. For instance:
        - Eigensolvers update `eigenvalues` and contribute to `kinetic`.
        - Potential solvers update `hartree`.
        - Exchange-correlation routines update `exchange` and `correlation` (and potentially `int_nvxc`).
        - Ion interaction calculations (e.g., in `esl_ion_interaction_m`) update `ionion`.
    - **Energy Consumers/Reporters:** After an SCF cycle or other energy calculation, the `energy_t` object holds the results. The `calculate` method is then called to sum these into `total`. The `display` method, or direct access to its members, is used for outputting results, logging, or for use in subsequent calculations (e.g., forces, MD steps).
    - **Physical Interpretation:** The specific formula in the `calculate` method is vital as it defines how the total energy is derived from its DFT components. The way various terms are added or subtracted (e.g., the subtraction of `hartree` and `int_nvxc`) is key to the correct physical interpretation of the total energy within the Kohn-Sham DFT framework.
