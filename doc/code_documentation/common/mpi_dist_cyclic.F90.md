## Overview

The `mpi_dist_cyclic_m` module implements a simple MPI-based cyclic data distribution scheme for one-dimensional arrays. It defines a derived type, `mpi_dist_cyclic_t`, which extends the generic `mpi_dist_t` base type (from `mpi_dist_m`). In this distribution, elements of a global array are assigned to MPI processes one by one in a round-robin fashion. For example, process 0 gets global element 1, process 1 gets global element 2, ..., process (P-1) gets global element P, process 0 gets global element P+1, and so on, where P is the number of processes.

This module provides the necessary functions to manage such a distribution, including initializing the scheme, determining the number of elements local to each process, and mapping between global and local indices, as well as identifying which process owns a particular global index. The functionality is primarily active when the code is compiled with MPI support (`WITH_MPI`).

*Note: The implementation details for some mapping functions in the source code appear to deviate from standard cyclic distribution algorithms, particularly concerning indexing and modulo operations. This documentation describes the implemented logic.*

## Key Components

- **Module:** `mpi_dist_cyclic_m`
    - **Description:** Implements a simple MPI-based cyclic distribution for 1D data arrays.
- **Type:** `mpi_dist_cyclic_t`, extends `mpi_dist_t`
    - **Description:** Represents a cyclic data distribution scheme. It inherits MPI communicator, rank, size, and global data size (`global_N`) from `mpi_dist_t`. It does not introduce new data components beyond those in the parent type.
    - **Procedures (Public):**
        - `init => init_`: Initializes the cyclic distribution.
        - `delete => delete_`: Resets the distribution object.
        - `local_N => local_N_`: Returns the number of elements local to the current process.
        - `glob_2_loc => global_2_local_`: Converts a global index to a local index (if owned).
        - `loc_2_glob => local_2_global_`: Converts a local index to its global index.
        - `glob_2_rank => global_2_rank_`: Determines the rank of the process owning a global index.

### Public Procedures (Type-Bound Methods or Functions)

#### `init_(this, comm, global_N)`
- **Description:** Initializes an `mpi_dist_cyclic_t` object. This routine calls the initializer of its parent type, `this%mpi_dist_t%init_(comm, global_N)`, which sets the MPI communicator, process rank, communicator size, and the total number of elements to be distributed. No actions specific to the cyclic nature (beyond what the parent `init_` does) are performed here.
- **Arguments:**
    - `this (class(mpi_dist_cyclic_t), intent(inout))`: The distribution object to initialize.
    - `comm (integer(im_), intent(in))`: The MPI communicator.
    - `global_N (integer(ii_), intent(in))`: The total number of elements in the global 1D array.

#### `delete_(this)`
- **Description:** Resets an `mpi_dist_cyclic_t` object by calling the `delete_` procedure of its parent type, `this%mpi_dist_t%delete_()`. This resets common MPI properties like `comm`, `rank`, `size`, and `global_N`.
- **Arguments:**
    - `this (class(mpi_dist_cyclic_t), intent(inout))`: The distribution object to reset.

#### `local_N_(this) result(N)`
- **Description:** Calculates the number of elements (`N`) that are stored locally on the calling MPI process. The implemented logic is: `N = this%global_N / this%size`. Then, it adjusts `N` by adding `ONE` (1) if `this%glob_2_rank(this%global_N) <= this%rank`. If compiled without MPI, it returns `this%global_N`.
*(Note: This logic for handling remainders is non-standard for typical cyclic distributions. A common approach is `N = global_N / size; if (rank < mod(global_N, size)) N = N + 1`.)*
- **Arguments:**
    - `this (class(mpi_dist_cyclic_t), intent(in))`: The distribution object.
- **Returns:** `N (integer(ii_))`: The number of elements local to the current process.

#### `global_2_local_(this, global_idx) result(local_idx)`
- **Description:** Converts a 1-based global index (`global_idx`) to its 1-based local index (`local_idx`) on the process that owns it. It first checks ownership using `this%glob_2_rank(global_idx)`. If the current process is the owner, the local index is computed as `global_idx / this%size` (integer division). If not owned, it returns `-ONE` (`-1_ii_`). If compiled without MPI, `local_idx = global_idx`.
*(Note: The local index calculation `global_idx / this%size` for 1-based indices is unusual. For a 0-based global index `g_idx`, the 0-based local index `l_idx` is typically `g_idx / number_of_processes`. For 1-based indices, this would translate to `(global_idx - 1) / number_of_processes + 1`.)*
- **Arguments:**
    - `this (class(mpi_dist_cyclic_t), intent(in))`: The distribution object.
    - `global_idx (integer(ii_), intent(in))`: The 1-based global index.
- **Returns:** `local_idx (integer(ii_))`: The 1-based local index if the element is local to the current process, otherwise `-1_ii_`.

#### `local_2_global_(this, local_idx) result(global_idx)`
- **Description:** Converts a 1-based local index (`local_idx`) on the calling MPI process to its corresponding 1-based global index (`global_idx`). The implemented formula is `global_idx = this%size * local_idx + this%rank`. If compiled without MPI, `global_idx = local_idx`.
*(Note: This formula appears to be for 0-based local and global indices. For 1-based local `l_idx` (from 1 to N_local) and 0-based rank `r`, a common 0-based global index `g_idx` is `r + (l_idx-1)*size`. The 1-based global index would then be `r + (l_idx-1)*size + 1`.)*
- **Arguments:**
    - `this (class(mpi_dist_cyclic_t), intent(in))`: The distribution object.
    - `local_idx (integer(ii_), intent(in))`: The 1-based local index on the current process.
- **Returns:** `global_idx (integer(ii_))`: The 1-based global index.

#### `global_2_rank_(this, global_idx) result(owner_rank)`
- **Description:** Determines the 0-based rank (`owner_rank`) of the MPI process that stores the element with the given 1-based global index (`global_idx`). The implemented formula is `owner_rank = mod(global_idx - ONE, this%global_N)`. If compiled without MPI, `owner_rank = 0`.
*(Note: This formula is incorrect for a standard cyclic distribution over `this%size` processes. It should be `owner_rank = mod(global_idx - ONE, this%size)`.)*
- **Arguments:**
    - `this (class(mpi_dist_cyclic_t), intent(in))`: The distribution object.
    - `global_idx (integer(ii_), intent(in))`: The 1-based global index.
- **Returns:** `owner_rank (integer(ii_))`: The 0-based rank of the MPI process.

## Important Variables/Constants
- **`ONE (integer(ii_), parameter)`**: A constant with value `1_ii_`, used in index calculations.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage:
!
! use mpi_dist_cyclic_m
! use prec, only: ii_
! #ifdef WITH_MPI
!   use mpi
! #endif
!
! type(mpi_dist_cyclic_t) :: cyclic_distribution
! integer :: my_rank, comm_size, global_data_size
! integer(ii_) :: num_local_elems, glob_idx, loc_idx, owner_process
! #ifdef WITH_MPI
!   integer :: comm_h = MPI_COMM_WORLD
!   call MPI_Comm_rank(comm_h, my_rank, ierr)
!   call MPI_Comm_size(comm_h, comm_size, ierr)
! #else
!   integer :: comm_h = 0
!   my_rank = 0
!   comm_size = 1
! #endif
!
! global_data_size = 17
! call cyclic_distribution%init(comm_h, int(global_data_size, kind=ii_))
!
! num_local_elems = cyclic_distribution%local_N()
! print *, "Rank ", my_rank, " has ", num_local_elems, " local elements (cyclic)."
!
! glob_idx = 7_ii_ ! Example global index
! owner_process = cyclic_distribution%global_2_rank(glob_idx)
! loc_idx = cyclic_distribution%global_2_loc(glob_idx) ! This is local index on 'owner_process'
!
! if (my_rank == owner_process) then
!    print *, "Rank ", my_rank, " owns global idx ", glob_idx, " at its local idx ", loc_idx
!    ! To verify:
!    ! print *, "  (Verification: local_2_global for this loc_idx: ", cyclic_distribution%local_2_global(loc_idx), ")"
! end if
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `mpi` (Conditional, if `WITH_MPI` is defined): Used by the parent type `mpi_dist_t` for MPI context (rank, size).
    - `mpi_dist_m`: Provides the parent type `mpi_dist_t`. This `mpi_dist_cyclic_t` inherits common MPI properties and initialization from it.
    - `prec, only: im_ => ip, ii_ => ip`: Imports integer kind parameters from the `prec` module.
- **External Libraries:**
    - **MPI Library** (Conditional, if `WITH_MPI` is defined): The underlying MPI library is essential for parallel execution and for the MPI context used by this distribution.
- **Interactions with other components:**
    - **Parent Type (`mpi_dist_t`):** `mpi_dist_cyclic_t` extends `mpi_dist_t` and uses its `init_` and `delete_` procedures for managing common MPI-related data.
    - **Parallel Data Management:** This type is intended for use by higher-level application components that need to distribute 1D arrays or data structures across MPI processes using a simple cyclic scheme.
    - **Algorithm Implementation:** Parallel algorithms that operate on data distributed in this manner would use the provided mapping functions (`local_N_`, `glob_2_loc`, `loc_2_glob`, `global_2_rank_`) to determine data locality, convert between local and global views, and manage communication.
    - **Build System:** The `WITH_MPI` preprocessor flag, controlled by the build system, dictates whether the parallel distribution logic is active or if the code falls back to a serial mode (where all data is considered local to rank 0).
    - **Correctness Concerns:** As noted, several of the implemented formulas for index mapping and local size calculation appear to deviate from standard, correct algorithms for 1-based cyclic distributions. Users should be cautious and verify these implementations if strict adherence to a standard cyclic distribution is required.
