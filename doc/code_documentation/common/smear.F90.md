## Overview

The `esl_smear_m` module is responsible for managing the calculation of electronic occupation numbers and the Fermi level, particularly for systems (like metals) where a smearing of the electronic states around the Fermi level is necessary for stable and efficient electronic structure calculations. It defines a derived type `smear_t` to hold smearing parameters (type of smearing, temperature/width) and the calculated Fermi level. The actual computation of occupations and the Fermi level is delegated to routines from the ELSI (Electronic Structure Library Interface) library, accessed via the `elsi` module.

The module supports several common smearing schemes, selectable through user input.

## Key Components

- **Module:** `esl_smear_m`
    - **Description:** Manages parameters and calculations related to Fermi-level smearing of electronic occupations.

- **Public Parameters (Smearing Scheme Identifiers):**
    These integer parameters are used to specify the type of smearing function:
    - `GAUSSIAN (integer, parameter, public)`: Value `0`. Represents Gaussian smearing.
    - `FD (integer, parameter, public)`: Value `1`. Represents Fermi-Dirac smearing.
    - `MP (integer, parameter, public)`: Value `2`. Represents Methfessel-Paxton smearing.
    - `CUBIC (integer, parameter, public)`: Value `3`. (Likely another type of smearing, e.g., cubic spline).
    - `COLD (integer, parameter, public)`: Value `4`. Represents Cold smearing by Marzari-Vanderbilt.

- **Type:** `smear_t`
    - **Description:** A derived type that encapsulates smearing parameters and the resulting Fermi level.
    - **Components:**
        - `smearing (integer(ip))`: An integer code (one of the parameters above) indicating the chosen smearing scheme.
        - `fermi_level (real(dp))`: The calculated Fermi level (chemical potential) in Hartree units.
        - `eTemp (real(dp))`: The smearing temperature in Hartree units. Default: `0.00095_dp` (approximately 300 Kelvin). Read from FDF key `Smearing.Temp`.
        - `eBroad (real(dp))`: The energy broadening width in Hartree units. Default: `0.1_dp`. Read from FDF key `eBroad`.
        - `occ_tol (real(dp))`: An occupation tolerance. (Declared in the source but not explicitly used in the provided subroutines; may be intended for future use or used internally by ELSI).
    - **Procedures (public, type-bound):**
        - `init`: Initializes the `smear_t` object by reading smearing scheme and parameters from an FDF input file.
        - `calc_fermi_occ`: Triggers the calculation of the Fermi level and updates electronic occupation numbers by interfacing with ELSI routines.

### Type-Bound Procedures

#### `init(this)`
- **Description:** Initializes the `smear_t` object `this`.
    1.  Reads a character string specifying the smearing scheme from the FDF input key `Smearing` (default is "FD" for Fermi-Dirac).
    2.  Converts this string identifier into the corresponding integer parameter (e.g., "FD" becomes `FD`) and stores it in `this%smearing`.
    3.  Reads the smearing temperature `this%eTemp` from the FDF key `Smearing.Temp` (defaults to `0.00095_dp` Hartree).
    4.  Reads the energy broadening width `this%eBroad` from the FDF key `eBroad` (defaults to `0.1_dp` Hartree).
- **Arguments:**
    - `this (class(smear_t), intent(inout))`: The `smear_t` object to be initialized.

#### `calc_fermi_occ(this, elsic, states)`
- **Description:** This subroutine orchestrates the calculation of the Fermi level and the electronic occupation numbers. It relies on external routines, typically from the ELSI library:
    1.  It calls `elsi_set_mu_broaden_scheme(elsic%e_h, this%smearing)` to pass the chosen smearing method (stored in `this%smearing`) to the ELSI library via the ELSI handle `elsic%e_h`.
    2.  It calls `elsi_set_mu_broaden_width(elsic%e_h, this%eBroad)` to pass the broadening width.
    3.  It then calls `elsi_compute_mu_and_occ(...)`. This ELSI routine takes the ELSI handle, total number of electrons (`states%nel`), details about the electronic states (number of states, spins, k-points, k-point weights, and the eigenvalues from `states%eigenvalues`), and computes:
        *   The occupation numbers, which are stored back into `states%occ_numbers`.
        *   The Fermi level, which is stored in `this%fermi_level`.
- **Arguments:**
    - `this (class(smear_t), intent(inout))`: The smearing object. Its `fermi_level` member is updated.
    - `elsic (type(elsi_t), intent(inout))`: An object of `elsi_t` (from `esl_elsi_m`), which provides the ELSI handle (`elsic%e_h`) needed for ELSI calls.
    - `states (type(states_t), intent(inout))`: An object of `states_t` (from `esl_states_m`). It provides input eigenvalues, k-point information, and total electron count to the ELSI routine, and its `occ_numbers` component is updated with the calculated occupations.

## Important Variables/Constants
- **Smearing Scheme Parameters:** `GAUSSIAN, FD, MP, CUBIC, COLD` are public integer parameters identifying the available smearing methods.
- **`smear_t` members:** `smearing`, `fermi_level`, `eTemp`, and `eBroad` are key components that define the smearing configuration and store the primary result (Fermi level).

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage within an SCF loop:
!
! module scf_cycle_m
!   use esl_smear_m, only: smear_t
!   use esl_elsi_m, only: elsi_t        ! For ELSI handle wrapper
!   use esl_states_m, only: states_t    ! For eigenvalues and occupations
!   implicit none
!
!   subroutine perform_scf_iteration(smear_config, elsi_handler, current_states)
!     type(smear_t), intent(inout) :: smear_config ! Initialized elsewhere
!     type(elsi_t), intent(inout) :: elsi_handler   ! Initialized elsewhere
!     type(states_t), intent(inout) :: current_states ! Eigenvalues are input, occupations are output
!
!     ! ... (Previous step: eigenvalues in current_states%eigenvalues are computed) ...
!
!     ! Calculate Fermi level and occupation numbers
!     call smear_config%calc_fermi_occ(elsi_handler, current_states)
!
!     ! Now, current_states%occ_numbers contains updated occupations,
!     ! and smear_config%fermi_level contains the new Fermi level.
!     ! These can be used to compute the new electron density and energy terms.
!
!   end subroutine perform_scf_iteration
! end module scf_cycle_m
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only: ip, dp`: Imports integer (`ip`) and double precision (`dp`) kind specifiers.
    - `fdf, only: fdf_get, leqi`: Used in the `init` procedure to read smearing scheme type and parameters from an FDF input file. `leqi` is used for case-insensitive string comparison.
    - `esl_states_m, only: states_t`: Provides the `states_t` derived type, which is essential for `calc_fermi_occ` as it supplies input eigenvalues and other state-related information, and receives the output occupation numbers.
    - `esl_elsi_m, only: elsi_t`: Provides the `elsi_t` type, which wraps the ELSI library handle (`e_h`) used in calls to ELSI routines.
    - `elsi` (module): This is the direct interface to the ELSI library (or its fake/stub version). `esl_smear_m` calls `elsi_set_mu_broaden_scheme`, `elsi_set_mu_broaden_width`, and `elsi_compute_mu_and_occ` from this `elsi` module.
- **External Libraries:**
    - **ELSI (Electronic Structure Library Interface):** The core functionality of calculating occupations and Fermi level is delegated to the ELSI library. This module acts as a high-level configuration and calling layer for ELSI's smearing-related capabilities.
- **Interactions with other components:**
    - **SCF Cycle (e.g., `esl_scf_m`):** The `calc_fermi_occ` routine is a critical part of the Self-Consistent Field (SCF) loop for metallic systems. After the electronic eigenvalues are computed in an iteration, `calc_fermi_occ` is called to determine the new Fermi level and update the occupation numbers of these states according to the chosen smearing scheme.
    - **Energy Calculation (e.g., `esl_energy_m`):** The calculated `fermi_level` is an important physical quantity often stored in an `energy_t` object. Furthermore, smearing schemes introduce an electronic entropy term to the free energy, which would be calculated using `this%eTemp` and the `states%occ_numbers`.
    - **Density Calculation (e.g., `esl_density_m`):** The updated `states%occ_numbers` are essential for computing the electron density for the next SCF iteration.
    - **Input Configuration (`fdf`):** The choice of smearing method and its associated parameters (temperature, broadening width) are made configurable by the user through an FDF input file.
    - **ELSI Interface (`elsi` module and `esl_elsi_m`):** This module directly uses the ELSI interface to perform its primary task. The `elsi_t` object from `esl_elsi_m` provides the necessary ELSI context (handle).
</tbody>
</table>
