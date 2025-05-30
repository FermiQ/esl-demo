## Overview

The `mpi_dist_m` module introduces a base derived type, `mpi_dist_t`, which is intended to serve as a foundational class for implementing various MPI (Message Passing Interface) based data distribution schemes. This base type is designed for distributing linear or one-dimensional quantities across multiple MPI processes. It stores essential MPI context, such as the communicator, the rank of the current process, the total number of processes in the communicator, and the global size of the data to be distributed.

The module itself does not implement specific distribution logic (like how global indices map to local indices or to specific ranks). Instead, it expects that derived types inheriting from `mpi_dist_t` will provide concrete implementations for functions like `glob_2_loc`, `loc_2_glob`, `glob_2_rank`, and `local_N`. The `mpi_dist_t` type provides basic, non-overridable initialization (`init_`) and reset (`delete_`) procedures.

## Key Components

- **Module:** `mpi_dist_m`
    - **Description:** Defines a generic base type `mpi_dist_t` for MPI-based distribution of one-dimensional data.
- **Type:** `mpi_dist_t`
    - **Description:** A base data type for MPI distributions. It encapsulates fundamental MPI properties and the global size of the distributed data.
    - **Components:**
        - `comm (integer(im_))`: The MPI communicator handle associated with this particular data distribution. Default value: `-1`. The kind `im_` is aliased from `ip` (likely a standard integer) from the `prec` module.
        - `rank (integer(ii_))`: The 0-based rank of the current MPI process within the `comm` communicator. Default value: `0`. The kind `ii_` is also aliased from `ip`.
        - `size (integer(ii_))`: The total number of MPI processes participating in the `comm` communicator. Default value: `1`.
        - `global_N (integer(ii_))`: The total number of elements in the one-dimensional data array that is being distributed across all processes. Default value: `0`.
    - **Procedures (public, non_overridable, type-bound):**
        - `init_`: Initializes the `mpi_dist_t` object by setting its components based on the provided communicator and global data size.
        - `delete_`: Resets the components of the `mpi_dist_t` object to their default, uninitialized states.

### Important Note on Derived Implementations
As stated in the source code comments, concrete distribution strategies are to be implemented in types that inherit from `mpi_dist_t`. These derived types are expected to provide procedures such as:
- `glob_2_loc`: Converts a global index to a local index for the current process.
- `loc_2_glob`: Converts a local index on the current process to its global index.
- `glob_2_rank`: Determines which MPI rank owns a given global index.
- `local_N`: Returns the number of elements stored locally on the current process.

### Type-Bound Procedures

#### `init_(this, comm, global_N)`
- **Description:** Initializes an `mpi_dist_t` object. It assigns the provided `comm` (MPI communicator) and `global_N` (total number of elements) to the corresponding members of `this`. If the code is compiled `WITH_MPI`:
    - It calls `MPI_Comm_rank(comm, this%rank, err)` to determine the rank of the current process.
    - It calls `MPI_Comm_rank(comm, this%size, err)` to determine the size of the communicator. *(Source code error: This should be `MPI_Comm_size(comm, this%size, err)`)*.
If compiled without `WITH_MPI`, `this%rank` is set to `0` and `this%size` is set to `1`.
- **Arguments:**
    - `this (class(mpi_dist_t), intent(inout))`: The `mpi_dist_t` object to be initialized.
    - `comm (integer(im_), intent(in))`: The MPI communicator to be used for this distribution.
    - `global_N (integer(ii_), intent(in))`: The total number of elements that this distribution will manage globally.

#### `delete_(this)`
- **Description:** Resets the members of the `mpi_dist_t` object to indicate an uninitialized or default state. It sets `this%comm = -1`, `this%rank = 0`, `this%size = 1`, and `this%global_N = 0`. This procedure does not deallocate or finalize the MPI communicator itself.
- **Arguments:**
    - `this (class(mpi_dist_t), intent(inout))`: The `mpi_dist_t` object to be reset.

## Important Variables/Constants
- The components of `mpi_dist_t` (`comm`, `rank`, `size`, `global_N`) are the core data managed by this base type.

## Usage Examples
```fortran
! TODO: Add usage example
! This is a base type. Usage would typically involve a derived type.
! Conceptual usage of a derived type 'my_dist_t' that extends 'mpi_dist_t':
!
! module my_specific_distribution_m
!   use mpi_dist_m, only: mpi_dist_t
!   use prec, only: ii_
! #ifdef WITH_MPI
!   use mpi
! #endif
!   implicit none
!
!   type, extends(mpi_dist_t) :: my_dist_t
!     ! ... specific components for this distribution ...
!   contains
!     procedure, public :: init
!     ! ... implementations for glob_2_loc, loc_2_glob, etc. ...
!   end type my_dist_t
!
! contains
!   subroutine init(self, mpi_communicator, total_elements)
!     class(my_dist_t), intent(inout) :: self
!     integer, intent(in) :: mpi_communicator ! Should be integer(im_)
!     integer(ii_), intent(in) :: total_elements
!
!     ! Call the parent constructor/initializer
!     call self%mpi_dist_t%init_(mpi_communicator, total_elements)
!
!     ! ... initialize specific parts of my_dist_t ...
!   end subroutine init
!
! end module my_specific_distribution_m
!
! ! Later, in some parallel region:
! ! use my_specific_distribution_m
! ! type(my_dist_t) :: distribution_scheme
! ! integer :: comm_world_handle, N_global_data
! ! #ifdef WITH_MPI
! !   comm_world_handle = MPI_COMM_WORLD
! ! #else
! !   comm_world_handle = 0 ! Or some other appropriate non-MPI value
! ! #endif
! ! N_global_data = 1000
! !
! ! call distribution_scheme%init(comm_world_handle, N_global_data)
! ! ! Now distribution_scheme can be used to manage data distribution
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `mpi` (Conditional, if `WITH_MPI` is defined): Provides the MPI API functions like `MPI_Comm_rank` and `MPI_Comm_size` (note the typo in the source for `size`).
    - `prec, only: im_ => ip, ii_ => ip`: Imports integer kind parameters from the `prec` module. `im_` is used for MPI communicator handles, and `ii_` is used for rank, size, and global counts. Both are aliased to `ip` (likely a standard integer kind) in this configuration.
- **External Libraries:**
    - **MPI Library** (Conditional, if `WITH_MPI` is defined): The underlying MPI library is required for parallel execution and to provide the functionality for the MPI calls.
- **Interactions with other components:**
    - **Derived Distribution Types:** This module is intended to be the parent class for more specific distribution types (e.g., `mpi_dist_block_cyclic_m`, `mpi_dist_cyclic_m`). These derived types will inherit from `mpi_dist_t`, call its `init_` procedure, and implement the actual logic for mapping global indices to local storage and MPI processes.
    - **Parallel Algorithms:** Any part of the application that deals with distributed arrays or parallel processing of 1D data would utilize objects of types derived from `mpi_dist_t` to understand data locality, perform index conversions, and determine communication patterns.
    - **MPI Environment:** The initialization logic directly interacts with the MPI environment (if `WITH_MPI` is active) to query the rank and size of the current process within a given communicator.
    - **Build System:** The `WITH_MPI` preprocessor flag, controlled by the build system, determines whether the MPI-specific code paths are enabled.
