## Overview

The `esl_ion_interaction_m` module is dedicated to computing the electrostatic interaction energy and the resulting forces between atomic ions (nuclei, typically treated as point charges) in a given system. It provides methods to handle both periodic systems, where long-range interactions necessitate techniques like Ewald summation, and isolated (non-periodic) systems, where a direct pairwise Coulomb summation is sufficient. The module defines a derived type `ion_interaction_t` to hold parameters relevant to these calculations, such as the Ewald summation parameter `alpha`.

## Key Components

- **Module:** `esl_ion_interaction_m`
    - **Description:** Calculates ion-ion electrostatic interaction energy and forces, with specific implementations for periodic (Ewald sum) and isolated (direct sum) systems.
- **Type:** `ion_interaction_t`
    - **Description:** A derived type to store parameters and provide methods for ion-ion interaction calculations.
    - **Components:**
        - `alpha (real(dp))`: The Ewald summation parameter, which balances the convergence of the real-space and reciprocal-space sums in the Ewald method. It is initialized to `0.21d0` by default but can be overridden via FDF input.
    - **Procedures (public, type-bound):**
        - `init`: Initializes the `ion_interaction_t` object, primarily setting the `alpha` parameter.
        - `calculate_periodic`: Computes ion-ion energy and forces for periodic systems using the Ewald summation technique.
        - `calculate_isolated`: Computes ion-ion energy and forces for isolated systems using direct Coulomb summation.

### Type-Bound Procedures

#### `init(this)`
- **Description:** Initializes the `ion_interaction_t` object. Its main action is to read the value for the Ewald summation parameter `this%alpha` from an FDF (Flexible Data Format) input file. It looks for the FDF key `Ewald.Alpha`. If the key is not found, `alpha` retains its default value of `0.21d0`.
- **Arguments:** `this (class(ion_interaction_t), intent(inout))`: The `ion_interaction_t` object to be initialized.

#### `calculate_periodic(this, geo, forces, eii)`
- **Description:** Calculates the ion-ion interaction energy (`eii`) and the forces on each atom (`forces`) for a system with periodic boundary conditions. This is achieved using the Ewald summation method, which splits the Coulomb sum into three parts:
    1.  **Real-space sum:** A rapidly converging sum over nearby ions and their periodic images, using the complementary error function (`erfc`) to screen interactions.
    2.  **Reciprocal-space sum:** A sum over reciprocal lattice vectors, accounting for the long-range part of the interaction. This involves structure factors `S(G) = sum_i Z_i * exp(i G . r_i)`.
    3.  **Self-interaction correction:** A term that subtracts the interaction of each ion with its own screening charge distribution.
    The forces are derived analytically from these energy terms.
- **Arguments:**
    - `this (class(ion_interaction_t), intent(in))`: Contains the Ewald parameter `alpha`.
    - `geo (type(geometry_t), intent(in))`: Provides all necessary geometric information: atomic charges (via `geo%species(idx)%z_ion`), atomic positions (`geo%xyz`), cell vectors (`geo%cell`), inverse cell matrix (`geo%icell`), and cell volume (`geo%vol`).
    - `forces (real(dp), intent(out) :: forces(:,:))`: A 2D array of shape `(3, natoms)` where the calculated force vector for each atom will be stored.
    - `eii (real(dp), intent(out))`: The calculated total ion-ion interaction energy.

#### `calculate_isolated(this, geo, forces, eii)`
- **Description:** Calculates the ion-ion interaction energy (`eii`) and forces (`forces`) for an isolated, non-periodic system. This is done by a straightforward pairwise summation of Coulomb interactions:
    - Energy: `eii = sum_{j>i} (Z_i * Z_j) / |r_i - r_j|`
    - Force on atom `i` from atom `j`: `F_ij = (Z_i * Z_j * (r_i - r_j)) / |r_i - r_j|^3`
    The total force on each atom is the sum of these pairwise forces.
- **Arguments:**
    - `this (class(ion_interaction_t), intent(in))`: Although `this` is an argument for interface consistency, its `alpha` member is not used in this routine.
    - `geo (type(geometry_t), intent(in))`: Provides atomic charges (via `geo%species(idx)%z_ion`) and positions (`geo%xyz`).
    - `forces (real(dp), intent(out) :: forces(:,:))`: A 2D array of shape `(3, natoms)` for storing the calculated forces.
    - `eii (real(dp), intent(out))`: The calculated total ion-ion interaction energy.

## Important Variables/Constants
- **`alpha (real(dp))`**: The Ewald summation parameter stored within the `ion_interaction_t` type. It is crucial for calculations in periodic systems.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage within a system update routine:
!
! subroutine update_system_ion_interactions(system_obj, is_periodic)
!   use esl_system_m, only: system_t
!   type(system_t), intent(inout) :: system_obj
!   logical, intent(in) :: is_periodic
!
!   if (is_periodic) then
!     call system_obj%ion_inter%calculate_periodic(system_obj%geo, &
!                                                 system_obj%force%ionion, &
!                                                 system_obj%energy%ionion)
!   else
!     call system_obj%ion_inter%calculate_isolated(system_obj%geo, &
!                                                  system_obj%force%ionion, &
!                                                  system_obj%energy%ionion)
!   end if
! end subroutine update_system_ion_interactions
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only : dp`: Imports the `dp` kind specifier for double precision real numbers.
    - `fdf, only: fdf_get`: Used in the `init` procedure to read the `Ewald.Alpha` parameter from an FDF input file.
    - `esl_constants_m`: Provides mathematical constants like `PI` (implicitly used via `sqrt(PI)` in `calculate_periodic`).
    - `esl_geometry_m`: Provides the `geometry_t` derived type, which supplies essential inputs like atomic charges (via species), positions, and cell parameters (for periodic calculations).
- **External Libraries:**
    - This module does not have direct dependencies on external libraries beyond standard Fortran intrinsics (like `erf`).
- **Interactions with other components:**
    - **`esl_system_m`:** An object of `ion_interaction_t` is typically a member of the main `system_t` object (e.g., `system%ion_inter`). The results of the calculations (energy `eii` and `forces`) are usually stored back into corresponding members of `system%energy` (e.g., `system%energy%ionion`) and `system%force` (e.g., `system%force%ionion`).
    - **`esl_geometry_m` (`geometry_t`):** This module is a heavy consumer of data from `geometry_t`, requiring atomic charges, positions, and, for periodic systems, cell vectors, inverse cell matrix, and cell volume.
    - **System Update Cycle (e.g., `system_t%update` or `esl_next_step_m`):** The `calculate_periodic` or `calculate_isolated` routines are typically called when the atomic geometry changes (e.g., during a molecular dynamics step or geometry optimization) to update the ion-ion contribution to the total energy and forces. These contributions are generally independent of the electronic structure (SCF cycle) itself, depending only on ionic positions and charges.
    - **`esl_force_m`:** The ion-ion forces calculated by this module are a specific component (`ionion`) of the total forces managed by the `force_t` type from `esl_force_m`.
    - **`esl_energy_m`:** The ion-ion energy calculated here contributes to the `ionion` field of the `energy_t` type from `esl_energy_m`.
