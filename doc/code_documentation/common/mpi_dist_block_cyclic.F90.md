## Overview

The `mpi_dist_block_cyclic_m` module implements a specific strategy for distributing one-dimensional data arrays across multiple MPI processes: the block-cyclic distribution. This method is widely used in parallel linear algebra libraries like ScaLAPACK. The module defines a derived type, `mpi_dist_block_cyclic_t`, which extends the generic `mpi_dist_t` base type (from `mpi_dist_m`).

In a block-cyclic distribution, data is divided into contiguous blocks of a specified size. These blocks are then distributed among the MPI processes in a cyclic (round-robin) manner. This module provides the necessary functions to manage such a distribution, including initializing the distribution scheme, determining the number of local elements on each process, and mapping between global and local indices, as well as identifying which process owns a particular global index. These functionalities are primarily active when the code is compiled with MPI support (`WITH_MPI`).

## Key Components

- **Module:** `mpi_dist_block_cyclic_m`
    - **Description:** Implements an MPI-based block-cyclic distribution for 1D data arrays, similar to ScaLAPACK.
- **Type:** `mpi_dist_block_cyclic_t`, extends `mpi_dist_t`
    - **Description:** Represents a block-cyclic data distribution scheme. It inherits MPI communicator, rank, size, and global data size (`global_N`) from `mpi_dist_t`.
    - **Components (in addition to inherited ones):**
        - `block (integer(ii_))`: The size of the contiguous blocks of data that are distributed cyclically among processes. Default value: `1_ii_`. (`ii_` is an integer kind from `prec`).
    - **Procedures (Public):**
        - `init => init_`: Initializes the block-cyclic distribution. (Overrides parent's concept).
        - `delete => delete_`: Resets the distribution object. (Overrides parent's concept).
        - `glob_2_loc => global_2_local_`: Converts a global index to a local index. (Implements expected interface).
        - `loc_2_glob => local_2_global_`: Converts a local index to its global index. (Implements expected interface).
        - `glob_2_rank => global_2_rank_`: Determines the rank of the process owning a global index. (Implements expected interface).
        - `local_N_` (effectively public, though not in `contains` with `=>`): A function that returns the number of elements local to the current process. (Implements expected interface `local_N`).

### Public Procedures (Type-Bound Methods or Functions)

#### `init_(this, comm, global_N, block_size)`
- **Description:** Initializes an `mpi_dist_block_cyclic_t` object. It first calls the initializer of its parent type, `this%mpi_dist_t%init_(comm, global_N)`, to set up the basic MPI context (communicator, rank, size) and the total number of elements. It then stores the specified `block_size` in `this%block`.
- **Arguments:**
    - `this (class(mpi_dist_block_cyclic_t), intent(inout))`: The distribution object to initialize.
    - `comm (integer(im_), intent(in))`: The MPI communicator over which the data is distributed.
    - `global_N (integer(ii_), intent(in))`: The total number of elements in the global 1D array.
    - `block_size (integer(ii_), intent(in))`: The size of the blocks to be distributed cyclically.

#### `delete_(this)`
- **Description:** Resets an `mpi_dist_block_cyclic_t` object to a default state. It sets `this%block` to `ONE` (which is `1_ii_`) and then calls the `delete_` procedure of its parent type, `this%mpi_dist_t%delete_()`, which resets `comm`, `rank`, `size`, and `global_N`.
- **Arguments:**
    - `this (class(mpi_dist_block_cyclic_t), intent(inout))`: The distribution object to reset.

#### `local_N_(this) result(N)`
- **Description:** Calculates the number of elements (`N`) that are stored locally on the calling MPI process under this block-cyclic distribution. The logic is adapted from ScaLAPACK's `numroc` (numeric routine for number of rows or columns). If the code is compiled without MPI (`#else` branch), it simply returns `this%global_N` (as all data is local).
- **Arguments:**
    - `this (class(mpi_dist_block_cyclic_t), intent(in))`: The distribution object.
- **Returns:** `N (integer(ii_))`: The number of elements local to the current process.

#### `global_2_local_(this, global_idx) result(local_idx)`
- **Description:** Converts a 1-based global index (`global_idx`) of an element in the distributed array to its 1-based local index (`local_idx`) on the process that owns it. The formula used is adapted from ScaLAPACK's `indxg2l` (index global to local). If compiled without MPI, it returns `local_idx = global_idx`.
- **Arguments:**
    - `this (class(mpi_dist_block_cyclic_t), intent(in))`: The distribution object.
    - `global_idx (integer(ii_), intent(in))`: The 1-based global index.
- **Returns:** `local_idx (integer(ii_))`: The 1-based local index.

#### `local_2_global_(this, local_idx) result(global_idx)`
- **Description:** Converts a 1-based local index (`local_idx`) on the calling MPI process to its corresponding 1-based global index (`global_idx`) in the overall distributed array. The formula is adapted from ScaLAPACK's `indxl2g` (index local to global). If compiled without MPI, it returns `global_idx = local_idx`.
- **Arguments:**
    - `this (class(mpi_dist_block_cyclic_t), intent(in))`: The distribution object.
    - `local_idx (integer(ii_), intent(in))`: The 1-based local index on the current process.
- **Returns:** `global_idx (integer(ii_))`: The 1-based global index.

#### `global_2_rank_(this, global_idx) result(owner_rank)`
- **Description:** Determines the 0-based rank (`owner_rank`) of the MPI process that stores the element with the given 1-based global index (`global_idx`). The formula is adapted from ScaLAPACK's `indxg2p` (index global to process). If compiled without MPI, it returns `owner_rank = 0`.
- **Arguments:**
    - `this (class(mpi_dist_block_cyclic_t), intent(in))`: The distribution object.
    - `global_idx (integer(ii_), intent(in))`: The 1-based global index.
- **Returns:** `owner_rank (integer(ii_))`: The 0-based rank of the MPI process that owns the element.

## Important Variables/Constants
- **`block (integer(ii_))`**: A component of `mpi_dist_block_cyclic_t` that defines the size of the blocks used in the cyclic distribution.
- **`ONE (integer(ii_), parameter)`**: A constant with value `1_ii_`, used in index calculations to handle 1-based vs. 0-based indexing conventions.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage:
!
! use mpi_dist_block_cyclic_m
! use prec, only: ii_
! #ifdef WITH_MPI
!   use mpi
! #endif
!
! type(mpi_dist_block_cyclic_t) :: distribution
! integer :: my_rank, comm_size, global_data_size, block_s
! integer(ii_) :: num_local_elements, i_global, i_local, owner_p
! #ifdef WITH_MPI
!   integer :: comm_handle = MPI_COMM_WORLD
!   call MPI_Comm_rank(comm_handle, my_rank, ierr)
!   call MPI_Comm_size(comm_handle, comm_size, ierr)
! #else
!   integer :: comm_handle = 0
!   my_rank = 0
!   comm_size = 1
! #endif
!
! global_data_size = 100
! block_s = 10
!
! call distribution%init(comm_handle, int(global_data_size, kind=ii_), int(block_s, kind=ii_))
!
! num_local_elements = distribution%local_N_()
! print *, "Rank ", my_rank, " has ", num_local_elements, " local elements."
!
! ! Example: Find owner of global index 25 and its local index there
! i_global = 25_ii_
! owner_p = distribution%global_2_rank(i_global)
! i_local = distribution%global_2_loc(i_global) ! This is local index on 'owner_p'
! if (my_rank == owner_p) then
!   print *, "Rank ", my_rank, " owns global index ", i_global, " at local index ", i_local
! end if
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `mpi` (Conditional, if `WITH_MPI` is defined): Required by the parent type `mpi_dist_t` for basic MPI operations (rank, size retrieval), although this module itself doesn't directly call MPI routines other than what the parent does. The logic heavily relies on MPI concepts.
    - `mpi_dist_m`: Provides the parent type `mpi_dist_t`, from which `mpi_dist_block_cyclic_t` inherits.
    - `prec, only: im_ => ip, ii_ => ip`: Imports integer kind parameters from the `prec` module. `im_` is used for MPI communicator handles (via parent), and `ii_` is used for ranks, sizes, global counts, and block sizes.
- **External Libraries:**
    - **MPI Library** (Conditional, if `WITH_MPI` is defined): The underlying MPI library is necessary for parallel execution and for the MPI context (communicator, rank, size) used by this distribution scheme.
- **Interactions with other components:**
    - **Parent Type (`mpi_dist_t`):** `mpi_dist_block_cyclic_t` extends `mpi_dist_t` and utilizes its `init_` and `delete_` procedures for managing common MPI properties.
    - **Parallel Data Structures:** This distribution type is intended to be used by higher-level modules or applications that need to manage one-dimensional arrays (or data that can be mapped to 1D) distributed across MPI processes in a block-cyclic manner. This is common for distributing vectors or the rows/columns of matrices in parallel numerical algorithms.
    - **ScaLAPACK-like Distribution:** The comments highlight that the distribution logic is adapted from ScaLAPACK. This suggests that data distributed using `mpi_dist_block_cyclic_t` would be laid out in memory in a way that is compatible with or similar to ScaLAPACK's 1D block-cyclic distribution, which is useful for interoperability or when implementing algorithms inspired by ScaLAPACK.
    - **Build System:** The `WITH_MPI` preprocessor flag, managed by the build system, determines whether the parallel distribution logic is active or if the code operates in a serial fallback mode (where all data is local to rank 0).
