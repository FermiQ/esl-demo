## Overview

The `esl_calc_density_matrix_ac_m` module serves as an interface to the ELSI (Electronic Structure Library Interface) library, specifically tailored for calculations involving atomic-centered (AC) basis sets where matrices (Hamiltonian, Overlap, Density Matrix) are represented in a sparse format. Its primary functions are to compute the density matrix (DM) and the energy-density matrix (EDM) by leveraging ELSI's capabilities. Additionally, it provides a routine to retrieve energy-related quantities like the Kohn-Sham energy, Fermi level, and electronic entropy from ELSI.

The functionality of this module is largely conditional on the `WITH_ELSI` preprocessor flag being defined during compilation. Some MPI-related setup is further conditional on `WITH_MPI`.

## Key Components

- **Module:** `esl_calc_density_matrix_ac_m`
    - **Description:** Interfaces with the ELSI library to compute the density matrix, energy-density matrix, and retrieve related energy terms, assuming sparse matrix representations for atomic-centered basis sets.

### Public Subroutines (Conditional on `WITH_ELSI`)

#### `set_elsi_sparsity_pattern_ac(elsi_ac, mpi_dist_bc_ac, sparsity_pattern_ac)`
- **Description:** (Available if `WITH_MPI` and `WITH_ELSI` are defined). Configures the ELSI library for handling distributed sparse matrices according to a block-cyclic distribution.
    1.  Sets the MPI communicator for ELSI using `elsi_set_mpi`, passing the communicator from `mpi_dist_bc_ac`.
    2.  Calculates the global total number of non-zero elements (`nnz_global`) by summing the local non-zero counts (`sparsity_pattern_ac%nz`) from all MPI processes using `MPI_Allreduce`.
    3.  Provides ELSI with the structure of the sparse matrices using `elsi_set_csc`. This includes `nnz_global`, local non-zero count, number of rows (`sparsity_pattern_ac%nr`), column indices (`sparsity_pattern_ac%column`), and row pointers (`sparsity_pattern_ac%rptr`), effectively describing the matrix in Compressed Sparse Column (CSC) or Row (CSR) format.
    4.  Sets the block size for the distributed sparse matrix representation in ELSI via `elsi_set_csc_blk`, using the block size from `mpi_dist_bc_ac`.
    5.  If the current process is rank 0, it sets the ELSI output verbosity level to 3 using `elsi_set_output` for debugging information.
- **Arguments:**
    - `elsi_ac (type(elsi_t), intent(inout))`: The ELSI wrapper object (from `esl_elsi_m`) containing the ELSI handle.
    - `mpi_dist_bc_ac (type(mpi_dist_block_cyclic_t), intent(in))`: Describes the block-cyclic distribution of data across MPI processes.
    - `sparsity_pattern_ac (type(sparse_pattern_t), intent(in))`: Defines the sparsity pattern of the matrices involved.
- **Returns:** Not applicable (Subroutine).

#### `calc_density_matrix_ac(elsi_ac, mat_H, mat_S, mat_DM)`
- **Description:** Computes the density matrix (`mat_DM`) by calling the ELSI routine `elsi_dm_real_sparse`. This routine typically solves the generalized eigenvalue problem or performs a similar calculation based on the provided sparse Hamiltonian matrix (`mat_H%M`) and sparse overlap matrix (`mat_S%M`). The resulting density matrix is stored in `mat_DM%M`. The Kohn-Sham energy computed by ELSI during this process is stored in `elsi_ac%KS_energy`. The current implementation assumes real arithmetic for the matrices.
- **Arguments:**
    - `elsi_ac (type(elsi_t), intent(inout))`: The ELSI wrapper object. Its `KS_energy` member is updated.
    - `mat_H (type(sparse_matrix_t), intent(inout))`: The sparse Hamiltonian matrix.
    - `mat_S (type(sparse_matrix_t), intent(inout))`: The sparse overlap matrix.
    - `mat_DM (type(sparse_matrix_t), intent(inout))`: The output sparse density matrix. Its non-zero values (`M` array) are populated by ELSI.
- **Returns:** Not applicable (Subroutine).

#### `calc_energy_density_matrix_ac(elsi_ac, mat_EDM)`
- **Description:** Computes the energy-weighted density matrix (`mat_EDM`) by calling the ELSI routine `elsi_get_edm_real_sparse`. This matrix is often used in calculating certain energy components, like the band energy. The result is stored in `mat_EDM%M`. The current implementation assumes real arithmetic.
- **Arguments:**
    - `elsi_ac (type(elsi_t), intent(inout))`: The ELSI wrapper object.
    - `mat_EDM (type(sparse_matrix_t), intent(inout))`: The output sparse energy-density matrix. Its non-zero values (`M` array) are populated by ELSI.
- **Returns:** Not applicable (Subroutine).

#### `get_energy_results(elsi_ac, KS_energy_out, fermi_level_out, entropy_out)`
- **Description:** Retrieves several key energy-related quantities that are typically computed and stored by ELSI during the solution process (e.g., density matrix calculation).
    1.  Calls `elsi_get_mu(elsi_ac%e_h, elsi_ac%fermi_level)` to fetch the Fermi level from ELSI and store it in the `fermi_level` member of the `elsi_ac` object.
    2.  Calls `elsi_get_entropy(elsi_ac%e_h, elsi_ac%entropy)` to fetch the electronic entropy and store it in the `entropy` member of `elsi_ac`.
    3.  The Kohn-Sham energy (`elsi_ac%KS_energy`, previously stored by `calc_density_matrix_ac`), the fetched Fermi level, and fetched entropy are then assigned to the output arguments `KS_energy_out`, `fermi_level_out`, and `entropy_out`, respectively.
- **Arguments:**
    - `elsi_ac (type(elsi_t), intent(inout))`: The ELSI wrapper object. Its `fermi_level` and `entropy` members are updated from ELSI.
    - `KS_energy_out (real(dp), intent(out))`: The Kohn-Sham energy.
    - `fermi_level_out (real(dp), intent(out))`: The Fermi level.
    - `entropy_out (real(dp), intent(out))`: The electronic entropy.
- **Returns:** Not applicable (Subroutine).

## Important Variables/Constants
- This module does not define any public Fortran `parameter` constants.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage within an SCF cycle for an AC basis calculation:
!
! module scf_ac_m
!   use esl_calc_density_matrix_ac_m
!   use esl_elsi_m, only: elsi_t
!   use esl_sparse_matrix_m, only: sparse_matrix_t
!   use esl_sparse_pattern_m, only: sparse_pattern_t
!   use mpi_dist_block_cyclic_m, only: mpi_dist_block_cyclic_t
!   implicit none
!
!   subroutine run_ac_scf_step(elsi_data, mpi_dist, sparsity, H_matrix, S_matrix, DM_matrix)
!     type(elsi_t), intent(inout) :: elsi_data
!     type(mpi_dist_block_cyclic_t), intent(in) :: mpi_dist
!     type(sparse_pattern_t), intent(in) :: sparsity
!     type(sparse_matrix_t), intent(inout) :: H_matrix, S_matrix, DM_matrix
!     real(dp) :: ks_e, mu_e, s_e
!
! #ifdef WITH_ELSI
! #  ifdef WITH_MPI
!     ! This would be called once at setup
!     ! call set_elsi_sparsity_pattern_ac(elsi_data, mpi_dist, sparsity)
! #  endif
!
!     ! In SCF iteration:
!     call calc_density_matrix_ac(elsi_data, H_matrix, S_matrix, DM_matrix)
!     call get_energy_results(elsi_data, ks_e, mu_e, s_e)
!     print *, "KS Energy: ", ks_e, " Fermi: ", mu_e, " Entropy: ", s_e
!
!     ! Optionally, calculate Energy Density Matrix
!     ! type(sparse_matrix_t) :: EDM_matrix
!     ! call EDM_matrix%init(sparsity) ! Assuming EDM has same sparsity
!     ! call calc_energy_density_matrix_ac(elsi_data, EDM_matrix)
!     ! ...
!     ! call EDM_matrix%delete()
! #endif
!   end subroutine run_ac_scf_step
! end module
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only: ip, dp`: Imports integer (`ip`) and double precision (`dp`) kind specifiers.
    - `esl_elsi_m, only: elsi_t`: Provides the `elsi_t` type, which is a wrapper around the ELSI library handle.
    - `esl_sparse_matrix_m, only: sparse_matrix_t`: Provides the `sparse_matrix_t` type for representing sparse matrices (H, S, DM, EDM).
    - `elsi` (module): This is the direct Fortran interface to the ELSI library (or its fake/stub version). The module uses several routines from it, such as `elsi_set_mpi`, `elsi_set_csc`, `elsi_dm_real_sparse`, etc.
    - `mpi_dist_block_cyclic_m, only: mpi_dist_block_cyclic_t` (Conditional on `WITH_MPI`): Used in `set_elsi_sparsity_pattern_ac` to describe the data distribution.
    - `esl_sparse_pattern_m, only: sparse_pattern_t` (Conditional on `WITH_MPI`): Used in `set_elsi_sparsity_pattern_ac` to define the matrix sparsity.
    - `mpi` (Conditional on `WITH_MPI`): Provides MPI routines like `MPI_Allreduce` and constants like `mpi_integer4`, `mpi_sum`.
- **External Libraries:**
    - **ELSI (Electronic Structure Library Interface):** This is the primary external library this module interfaces with. All core computations (DM, EDM, energies) are delegated to ELSI.
    - **MPI (Message Passing Interface):** Required if the `WITH_MPI` flag is set, for parallel setup and execution.
- **Interactions with other components:**
    - **Hamiltonian and Overlap Construction (e.g., `esl_hamiltonian_ac_m`, `esl_overlap_matrix_ac_m`):** These modules are responsible for building the sparse Hamiltonian (`mat_H`) and overlap (`mat_S`) matrices that serve as input to `calc_density_matrix_ac`.
    - **SCF Cycle (e.g., `esl_scf_m` for AC basis):** This module is a key part of the SCF cycle when using ELSI for solving the eigensystem (or directly obtaining the DM).
        - `set_elsi_sparsity_pattern_ac` is typically called once during initialization if MPI is used.
        - `calc_density_matrix_ac` is called in each SCF iteration to get the updated density matrix.
        - `get_energy_results` is called to retrieve energy terms and the Fermi level, which are used for checking convergence and calculating the total energy.
    - **Density Management (e.g., `esl_density_matrix_ac_m`):** The density matrix (`mat_DM`) computed here is a central quantity. It would be managed by, or its data transferred to, a main density object used by other parts of the code (e.g., to calculate the Hartree and XC potentials for the next SCF iteration).
    - **Build System:** The `WITH_ELSI` and `WITH_MPI` preprocessor flags, managed by the build system, determine whether this module's functionalities are compiled and available.
</tbody>
</table>
