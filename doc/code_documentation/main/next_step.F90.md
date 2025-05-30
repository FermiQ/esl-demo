## Overview

The `esl_next_step_m` module provides a crucial subroutine, `next_step_setup`, which is designed to be invoked at the beginning of a new major computational segment, typically one that precedes a Self-Consistent Field (SCF) calculation. This module acts as a preparatory stage for various simulation events such as a Molecular Dynamics (MD) step, a step in a geometry optimization, or adjustments to fundamental computational parameters like energy cutoffs. Its current implementation includes a hook for Lua scripting and calls to update the system state based on the chosen basis set (Plane-Wave or Atomic-Centered).

## Key Components

- **Module:** `esl_next_step_m`
    - **Description:** This module's primary purpose is to manage the setup and preliminary operations required before embarking on a new significant computational step, which usually revolves around an SCF cycle.
    - **Public Procedures:**
        - `next_step_setup`

### Subroutine/Function: `next_step_setup(system)`

- **Description:** This subroutine orchestrates the necessary preparations for the subsequent computational step. Its actions include:
    1.  Invoking a Lua script hook identified by `LUA_INIT_STEP`, if the `WITH_FLOOK` preprocessor directive is enabled. This allows for external script-based modifications or initializations at the beginning of a step.
    2.  Branching its execution based on the type of basis set specified in `system%basis%type`:
        *   If Plane-Waves (`PLANEWAVES`): It calls the internal subroutine `next_planewave()`.
        *   If Atomic-Centered orbitals (`ATOMCENTERED`): It calls the internal subroutine `next_atomicorbs()`.
- **Arguments:**
    - `system (type(system_t), intent(inout))`: The main system object (`system_t` from `esl_system_m`). This object is passed with `intent(inout)` suggesting it can be modified by this setup routine, for instance, by updating atomic positions in an MD step.
- **Internal Subroutines:**
    - `next_atomicorbs()`: This subroutine is meant to contain specific setup logic for when an atomic-centered basis is used. Currently, its main action is to call `system%update(periodic=.false.)`.
    - `next_planewave()`: This subroutine is intended for setup logic specific to plane-wave basis calculations. Similar to `next_atomicorbs`, it currently calls `system%update(periodic=.false.)`.

## Important Variables/Constants

- `LUA` (Module variable from `esl_flook_global_m`, used if `WITH_FLOOK` is defined): This variable serves as a handle or state representation for the Lua interpreter, enabling interaction with Lua scripts.

## Usage Examples

```fortran
! TODO: Add usage example
! This subroutine is typically called from a main loop before starting an SCF cycle.
! Example (conceptual, within a larger loop in a main program):
!
! use esl_next_step_m
! use esl_system_m
! type(system_t) :: my_system
!
! ! ... (initialize my_system) ...
!
! do istep = 1, max_steps
!   call next_step_setup(my_system) ! Prepare for the new step (e.g., update geometry)
!
!   ! ... (perform SCF calculation using my_system) ...
! end do
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `esl_flook_global_m`, `esl_flook_if_m` (Conditional, if `WITH_FLOOK` is defined): These modules provide the necessary infrastructure for Lua scripting, including the `LUA` variable and the `flook_if_call` function to execute Lua scripts at the `LUA_INIT_STEP` hook.
    - `esl_system_m, only: system_t`: This provides the definition of the `system_t` derived type, which encapsulates the overall state of the simulated system. The `system%update()` method is a key interaction point.
    - `esl_basis_m, only: PLANEWAVES, ATOMCENTERED`: This module supplies the constants `PLANEWAVES` and `ATOMCENTERED`, used to determine the type of basis set and direct the flow within `next_step_setup`.
- **External Libraries:**
    - `Lua` (Conditional, if `WITH_FLOOK` is defined): The Lua scripting language library would be linked if this feature is enabled.
- **Interactions with other components:**
    - **`esl_main` (or similar driver routine):** `next_step_setup` is typically invoked by a higher-level control routine (like `esl_main` in this project) at the beginning of each major iteration or step (e.g., before each SCF cycle in an MD simulation).
    - **`esl_system_m`:** The most direct interaction is with the `system_t` object from `esl_system_m`. The call to `system%update()` suggests that `next_step_setup` is responsible for triggering updates to the system's state (e.g., atomic positions, cell parameters) that are necessary before the next SCF calculation can proceed.
    - **Lua Scripts (if `WITH_FLOOK` is enabled):** The `LUA_INIT_STEP` hook allows external Lua scripts to interface with the simulation, potentially modifying system parameters, logging data, or controlling the simulation flow at the start of a new step.
    - The comments within the module ("`Move quantities required for a move`", "`MD step`") strongly indicate that a primary intended use case for this module is in the context of Molecular Dynamics simulations, where atomic positions are updated at each step, requiring subsequent recalculation of electronic structure.
