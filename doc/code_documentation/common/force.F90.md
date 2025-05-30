## Overview

The `esl_force_m` module introduces a derived type, `force_t`, specifically designed for storing and managing atomic forces within an electronic structure simulation. Atomic forces are crucial for applications like geometry optimization and molecular dynamics. This module provides a structure to hold different contributions to the total force acting on each atom, such as local potential contributions, non-local (e.g., pseudopotential) contributions, and direct ion-ion interaction forces. All forces are expressed in Hartree/Bohr atomic units. The module includes routines for initializing these force arrays, calculating the total force from its components, displaying the force information, and ensuring proper memory deallocation when a `force_t` object is destroyed.

## Key Components

- **Module:** `esl_force_m`
    - **Description:** Defines the `force_t` data type for managing various components of atomic forces and provides associated procedures.
- **Type:** `force_t`
    - **Description:** A derived type that serves as a container for different force contributions acting on the atoms in the system. Each force component is represented as a 2D allocatable array of `real(dp)`, with dimensions `(3, natoms)`, where the first index corresponds to the Cartesian components (x, y, z) and the second index identifies the atom.
    - **Components (all `real(dp), allocatable :: name(:,:)`):**
        - `total`: Stores the total force vector for each atom, calculated as the sum of all other defined contributions.
        - `loc`: Stores the force contributions arising from local potentials (e.g., the local part of pseudopotentials, external local fields).
        - `nl`: Stores the force contributions from non-local potentials, most commonly the non-local part of pseudopotentials (like Kleinman-Bylander projectors).
        - `ionion`: Stores the forces due to direct electrostatic interactions between the atomic ions (nuclei).
    - **Procedures (public, type-bound):**
        - `init`: Allocates and initializes the force arrays.
        - `calculate`: Computes the `total` force by summing the `loc`, `nl`, and `ionion` components.
    - **Finalizer (type-bound):**
        - `finalizer`: Automatically deallocates the memory of the force arrays when a `force_t` object is finalized.

### Module Subroutine (Public)
#### `display(this, only_total)`
- **Description:** Prints the force components stored in a `force_t` object using YAML formatting.
- **Arguments:**
    - `this (class(force_t), intent(in))`: The `force_t` object whose data is to be displayed.
    - `only_total (logical, intent(in), optional)`: If present and true, only the `total` force array is printed. If absent or false, a detailed breakdown including `total`, `loc`, `nl`, and `ionion` forces is printed. *(Note: The source code's default behavior if `only_total` is absent is to print only the total force. For a full breakdown, `only_total` must be explicitly passed as `.false.`)*.

### Type-Bound Procedures

#### `init(this, natoms)`
- **Description:** Initializes a `force_t` object. This involves allocating memory for each of the force arrays (`total`, `loc`, `nl`, `ionion`) with dimensions `(3, natoms)`. All elements of these newly allocated arrays are set to `0.0_dp`.
- **Arguments:**
    - `this (class(force_t), intent(inout))`: The `force_t` object to be initialized.
    - `natoms (integer, intent(in))`: The number of atoms in the system, which dictates the second dimension of the force arrays.

#### `calculate(this)`
- **Description:** Computes the `this%total` force array by performing an element-wise sum of the `this%loc`, `this%nl`, and `this%ionion` arrays. This routine should be called after all individual force contributions have been calculated and stored in their respective arrays.
- **Arguments:**
    - `this (class(force_t), intent(inout))`: The `force_t` object. Its `total` component will be updated based on the sum of `loc`, `nl`, and `ionion`.

#### `finalizer(this)`
- **Description:** This is a type-bound `final` subroutine, automatically invoked by the Fortran runtime when an object of type `force_t` is deallocated or goes out of scope. It checks if the `total`, `loc`, `nl`, and `ionion` arrays are allocated and, if so, deallocates them to prevent memory leaks.
- **Arguments:**
    - `this (type(force_t), intent(inout))`: The `force_t` object being finalized.

## Important Variables/Constants
- The array components of `force_t` (`total`, `loc`, `nl`, `ionion`) are the primary entities, storing Cartesian force vectors for each atom.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage:
!
! use esl_force_m
! use esl_geometry_m ! Assuming geometry_t contains n_atoms
! type(force_t) :: atomic_forces
! type(geometry_t) :: sys_geometry
! integer :: num_atoms
!
! ! ... (initialize sys_geometry and get num_atoms) ...
! num_atoms = sys_geometry%n_atoms
!
! ! Initialize force object
! call atomic_forces%init(num_atoms)
!
! ! ... (During/after an SCF calculation, populate individual force components) ...
! ! e.g., atomic_forces%loc(:,atom_idx) = calculated_local_force_on_atom_idx
! !       atomic_forces%nl(:,atom_idx) = calculated_nonlocal_force_on_atom_idx
! !       atomic_forces%ionion(:,atom_idx) = calculated_ionion_force_on_atom_idx
!
! ! Calculate the total forces
! call atomic_forces%calculate()
!
! ! Display forces (e.g., only total force)
! call display(atomic_forces, only_total=.true.)
!
! ! Display all force components
! call display(atomic_forces, only_total=.false.)
!
! ! atomic_forces will be automatically finalized when it goes out of scope.
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only : dp`: Imports the `dp` kind specifier for double precision real numbers, used for all force arrays.
    - `yaml_output` (used by the module subroutine `display`): This module is required for generating YAML-formatted output of the force data.
- **External Libraries:**
    - This module has no direct dependencies on external libraries.
- **Interactions with other components:**
    - **`esl_system_m`:** Typically, an instance of `force_t` is a member of a higher-level system object (e.g., `system_t` from `esl_system_m`, accessed like `system%force`). This makes `force_t` an integral part of the overall system state.
    - **Force Calculation Engines:** Various specialized routines within the electronic structure code are responsible for computing the different force contributions. For example:
        - Hellmann-Feynman theorem components would contribute to `loc` and `nl`.
        - Pulay forces (if applicable) would also be added to these.
        - Ion-ion interaction routines (like those in `esl_ion_interaction_m`) would populate `ionion`.
    - **SCF Solvers / DFT Modules:** These modules usually orchestrate the calculation of electronic structure and then trigger the computation of forces. The resulting forces are stored in a `force_t` object.
    - **Geometry Optimizers and Molecular Dynamics Engines:** The `total` force array computed by `atomic_forces%calculate()` is the primary input for algorithms that update atomic positions, such as those found in geometry optimization schemes or molecular dynamics simulations.
    - **Output and Logging:** The `display` subroutine provides a standard way to report calculated forces, which is useful for debugging, analysis, and user feedback.
    - **Memory Management:** The `finalizer` subroutine ensures that the dynamically allocated memory for the force arrays is correctly released, which is crucial for long or multiple calculations to prevent memory leaks.
