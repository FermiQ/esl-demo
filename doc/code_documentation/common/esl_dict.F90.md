## Overview

The `esl_dict_m` module, when compiled with the `WITH_FLOOK` preprocessor flag, provides a mechanism for creating and populating a global dictionary named `esl_variables`. This dictionary serves as a bridge, allowing Fortran program variables to be exposed to an external scripting environment, presumably Lua through the Flook library. The module defines a generic interface `esl_dict_var_add` which simplifies the process of adding variables of different data types (character strings, logicals, integers, double precision reals) and ranks (scalars, 1D arrays, 2D arrays for reals) into this shared dictionary.

If `WITH_FLOOK` is not defined, this module has no active code.

## Key Components

- **Module:** `esl_dict_m`
    - **Description:** Implements a global dictionary (`esl_variables`) for sharing Fortran variables with a scripting environment (Flook/Lua). This functionality is conditional on the `WITH_FLOOK` preprocessor flag.
- **Global Variable (Conditional on `WITH_FLOOK`):**
    - `esl_variables (type(dict), public)`: A public dictionary object provided by the Flook `dictionary` module. It stores key-value pairs, where keys are character strings (variable names) and values can be Fortran scalars or arrays of various types.
- **Generic Interface (Conditional on `WITH_FLOOK`):**
    - `esl_dict_var_add`: A public generic interface that provides a unified way to add different types of Fortran variables to the `esl_variables` dictionary. It resolves to specific internal subroutines based on the data type and rank of the variable being added.
    - **Specific Module Procedures for `esl_dict_var_add`:**
        - `dict_variable_add_v_0d(name, val)`: Adds a character string `val`. The string is trimmed and stored by value (`.kv.`). If `name` already exists, the old string value in the dictionary is deallocated.
        - `dict_variable_add_b_0d(name, val)`: Adds a scalar logical `val`. The variable is passed with `intent(inout), target` and added as a pointer (`.kvp.`).
        - `dict_variable_add_i_0d(name, val)`: Adds a scalar integer `val`. Passed with `intent(inout), target` and added as a pointer (`.kvp.`).
        - `dict_variable_add_i_1d(name, val)`: Adds a 1D integer array `val`. Passed with `intent(inout), target` and added as a pointer (`.kvp.`).
        - `dict_variable_add_d_0d(name, val)`: Adds a scalar `real(dp)` `val`. Passed with `intent(inout), target` and added as a pointer (`.kvp.`).
        - `dict_variable_add_d_1d(name, val)`: Adds a 1D `real(dp)` array `val`. Passed with `intent(inout), target` and added as a pointer (`.kvp.`).
        - `dict_variable_add_d_2d(name, val)`: Adds a 2D `real(dp)` array `val`. Passed with `intent(inout), target` and added as a pointer (`.kvp.`).
    - **Update Behavior:** For all procedures, if a variable with the given `name` already exists in `esl_variables`, the old entry is removed before adding the new one. For the string version (`_v_0d`), the old dictionary string is deallocated. For pointer versions (`.kvp.`), only the dictionary's pointer entry is removed; the Fortran variable itself is not deallocated by this module.

## Important Variables/Constants
- **`esl_variables (type(dict))`**: The global dictionary instance itself. This is the central object managed by this module.

## Usage Examples
```fortran
! TODO: Add usage example
! This module is typically used by other Fortran modules to expose variables to Lua.
! Example from another module (e.g., esl_system_m%init):
!
! #ifdef WITH_FLOOK
!   use esl_dict_m, only: esl_dict_var_add
!   ! ... other declarations ...
!   real(dp) :: total_energy_val
!   integer :: num_atoms_val
!   character(len=20) :: system_name_val
!   real(dp), dimension(3, num_atoms_val), target :: atomic_coords_val
!
!   ! ... (assign values to total_energy_val, num_atoms_val, system_name_val, atomic_coords_val) ...
!
!   call esl_dict_var_add('E.Total', total_energy_val)
!   call esl_dict_var_add('Geometry.NumAtoms', num_atoms_val)
!   call esl_dict_var_add('System.Name', system_name_val)
!   call esl_dict_var_add('Geometry.xyz', atomic_coords_val) ! Assuming atomic_coords_val is suitably dimensioned
! #endif
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only : dp`: Imports the `dp` kind specifier for double precision real numbers.
    - `dictionary` (Conditional, if `WITH_FLOOK` is defined): This is a core dependency from the Flook library. It provides the `dict` type and the operators/methods used for dictionary manipulations (e.g., `.in.`, `delete`, `//`, `.kv.` for key-value copy, `.kvp.` for key-value pointer).
- **External Libraries:**
    - **Flook/Lua** (Implicit, if `WITH_FLOOK` is defined): The functionality of this module is intrinsically linked to the Flook library, which facilitates Fortran-Lua interoperability. The `dictionary` module is a component of Flook.
- **Interactions with other components:**
    - **Provider of Variables to Lua:** Various Fortran modules within the application (e.g., `esl_system_m`, `esl_scf_m`) call `esl_dict_var_add` to register their key variables. These variables then become accessible from the Lua scripting environment via Flook.
    - **Scripting Environment (`esl_flook_if_m`, `esl_flook_global_m`):** Other Flook-related modules like `esl_flook_if_m` or `esl_flook_global_m` would be responsible for initializing the Lua interpreter and providing the Lua-side mechanisms to interact with the `esl_variables` dictionary.
    - **Data Synchronization (for pointer variables):** When variables are added using `.kvp.` (which applies to logicals, integers, and reals and their arrays), the dictionary stores a pointer to the original Fortran variable. This means:
        - If Lua modifies the variable via the dictionary, the Fortran variable's value changes.
        - If the Fortran code modifies the variable after it has been added to the dictionary, Lua scripts accessing it will see the updated value.
    - **Data Copying (for string variables):** For character strings added with `dict_variable_add_v_0d` (using `.kv.`), a copy of the string's value is stored in the dictionary. Subsequent changes to the original Fortran string variable will not affect the value in the dictionary, and vice-versa.
- **Conditional Compilation:** The entire functionality of this module is contingent on the `WITH_FLOOK` preprocessor flag. If this flag is not set during compilation, the module will be effectively empty, and no dictionary operations will be available.
