## Overview

The `elsi` module provided in the file `elsi_fake.F90` serves as a **stub** or **mock implementation** of the ELSI (Electronics Structure Library Interface) API. Its primary purpose is to allow the main electronic structure application, which is designed to interface with the ELSI library, to compile and run even when the actual ELSI library is not available or not linked. This is extremely useful for development, basic testing, or for running on systems where ELSI is not installed.

While many of the ELSI functionalities related to solving eigenvalue problems or performing complex matrix operations are implemented as empty shells (i.e., they do nothing), this fake module provides a working implementation for calculating the Fermi level and electronic occupation numbers (`elsi_compute_mu_and_occ`) using a Fermi-Dirac smearing scheme.

## Key Components

- **Module:** `elsi`
    - **Description:** A substitute implementation for the ELSI library interface, providing placeholder routines and a functional Fermi-Dirac smearing calculation.
- **Type:** `elsi_handle`
    - **Description:** A minimal derived type designed to mimic the opaque handle used by the actual ELSI library. In this fake implementation, it's used to store parameters relevant to the Fermi level calculation, such as the number of electrons, the chosen smearing method, and the smearing width.
    - **Components:**
        - `nel (real(dp))`: Stores the total number of electrons.
        - `method (integer(ip))`: An integer code representing the selected smearing/broadening scheme (e.g., `1` for Fermi-Dirac).
        - `width (real(dp))`: The energy broadening width used in the smearing function.

### Subroutine/Function: `elsi_init(eh, solver, mode, mat_format, nbasis, nel, nstate)`
- **Description:** A fake implementation of the ELSI initialization routine. In this version, it only captures the number of electrons (`nel`) and stores it in the `elsi_handle` (`eh`). All other parameters (like `solver`, `mode`, `mat_format`, `nbasis`, `nstate`) are ignored.
- **Arguments:**
    - `eh (type(elsi_handle), intent(inout))`: The ELSI handle to be "initialized."
    - `nel (real(dp), intent(in))`: The total number of electrons.
    - `solver, mode, mat_format, nbasis, nstate (various types, intent(in))`: These arguments are accepted but not used by the fake routine.

### Subroutine/Function: `elsi_finalize(eh)`
- **Description:** A fake implementation of the ELSI finalization routine. It performs no actions.
- **Arguments:** `eh (type(elsi_handle), intent(inout))`: The ELSI handle.

### Subroutine/Function: `elsi_set_csc(eh, global_nnz, local_nnz, nrow, col_idx, row_ptr)`
- **Description:** A fake routine that would normally pass Compressed Sparse Column (CSC) matrix data to ELSI. This stub version performs no actions.
- **Arguments:** All arguments are accepted but ignored.

### Subroutine/Function: `elsi_dm_real_sparse(eh, h, s, d, energy)`
- **Description:** A fake routine that would typically calculate the density matrix (`d`) from Hamiltonian (`h`) and overlap (`s`) matrices for real sparse systems using ELSI. This stub performs no calculation and does not modify `d` or `energy`.
- **Arguments:** `h, s, d (real(dp), intent(inout))`, `energy (real(dp), intent(out))`.

### Subroutine/Function: `elsi_get_edm_real_sparse(eh, edm)`
- **Description:** A fake routine that would get the energy-density matrix (`edm`) from ELSI. This stub does not modify `edm`.
- **Arguments:** `edm (real(dp), intent(inout))`.

### Subroutine/Function: `elsi_compute_mu_and_occ(eh, nel, nstate, nspin, nkpt, k_weights, evals, occ, mu)`
- **Description:** This is the most significant functional part of the fake ELSI module. It calculates the Fermi level (`mu`) and electronic occupation numbers (`occ`) based on input eigenvalues (`evals`), the total number of electrons (`nel`), and smearing parameters (method and width) stored in the `elsi_handle` (`eh`). The calculation of `mu` is performed using a bisection method to find the chemical potential that yields the correct total number of electrons when summed over states using the specified smearing function (currently, only Fermi-Dirac is implemented).
- **Internal Subroutines:**
    - `next_step(mu, diff)`: A helper that selects the smearing function based on `eh%method`. Currently, it only proceeds if `eh%method` is 1 (Fermi-Dirac).
    - `fermi_dirac(mu, diff)`: Calculates occupation numbers for all states using the Fermi-Dirac distribution for a given chemical potential `mu` and broadening `eh%width`. It also computes `diff`, the difference between the sum of these occupations (weighted by `k_weights`) and the target number of electrons `nel`.
- **Arguments:**
    - `eh (type(elsi_handle), intent(in))`: The ELSI handle containing smearing parameters.
    - `nel (real(dp), intent(in))`: Target total number of electrons.
    - `nstate, nspin, nkpt (integer(ip), intent(in))`: Dimensions: number of states, spins, and k-points.
    - `k_weights (real(dp), dimension(nkpt), intent(in))`: Weights for each k-point.
    - `evals (real(dp), dimension(nstate,nspin,nkpt), intent(in))`: Eigenvalues for each state, spin, and k-point.
    - `occ (real(dp), dimension(nstate,nspin,nkpt), intent(out))`: Output array for calculated occupation numbers.
    - `mu (real(dp), intent(out))`: Output calculated Fermi level.

### Subroutine/Function: `elsi_set_mu_broaden_scheme(eh, method)`
- **Description:** Sets the smearing/broadening scheme to be used by `elsi_compute_mu_and_occ`. The integer `method` is stored in `eh%method`.
- **Arguments:**
    - `eh (type(elsi_handle), intent(inout))`: The ELSI handle.
    - `method (integer(ip), intent(in))`: An integer identifying the broadening scheme (e.g., 1 for Fermi-Dirac).

### Subroutine/Function: `elsi_set_mu_broaden_width(eh, width)`
- **Description:** Sets the energy broadening width for the smearing scheme. The `width` is stored in `eh%width`.
- **Arguments:**
    - `eh (type(elsi_handle), intent(inout))`: The ELSI handle.
    - `width (real(dp), intent(in))`: The energy broadening width.

## Important Variables/Constants
- This module does not define any public Fortran `parameter` constants. The `elsi_handle` type internally stores `nel`, `method`, and `width` for the stub operations.

## Usage Examples
```fortran
! TODO: Add usage example
! This module is used implicitly when the main code calls ELSI routines.
! If this fake module is linked instead of the real ELSI library, its routines are called.

! Example of setting up smearing and getting Fermi level:
! use elsi
! use prec, only: dp, ip
! type(elsi_handle) :: my_elsi_handle
! real(dp) :: num_electrons, broadening_width, fermi_energy
! real(dp), allocatable :: eigenvalues(:,:,:), occupations(:,:,:)
! real(dp), allocatable :: kpt_weights(:)
! integer :: n_states, n_spins, n_kpts, smearing_method_id
!
! ! ... (initialize n_states, n_spins, n_kpts, num_electrons, broadening_width) ...
! ! ... (allocate and fill eigenvalues, kpt_weights) ...
! ! ... (allocate occupations) ...
!
! smearing_method_id = 1 ! Fermi-Dirac
! call elsi_init(my_elsi_handle, 0, 0, 0, 0, num_electrons, 0) ! Dummy values for non-nel args
! call elsi_set_mu_broaden_scheme(my_elsi_handle, smearing_method_id)
! call elsi_set_mu_broaden_width(my_elsi_handle, broadening_width)
!
! call elsi_compute_mu_and_occ(my_elsi_handle, num_electrons, &
!                             n_states, n_spins, n_kpts, kpt_weights, &
!                             eigenvalues, occupations, fermi_energy)
!
! ! ... (use occupations and fermi_energy) ...
!
! call elsi_finalize(my_elsi_handle)
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - `prec, only: ip, dp`: Imports integer (`ip`) and double precision (`dp`) kind specifiers from the `prec` module. These are used for variable declarations.
- **External Libraries:**
    - None. This module specifically avoids the need for the actual external ELSI library by providing these substitute routines.
- **Interactions with other components:**
    - **Transparent Replacement:** This module is designed to be a drop-in replacement for the actual ELSI interface. Other modules in the application (e.g., `esl_smear_m` which calculates occupations, or solver routines in `esl_hamiltonian_m` or `esl_scf_m` which would call ELSI for matrix solutions) would `use elsi`. If this `elsi_fake.F90` is compiled and linked, its versions of the subroutines are called.
    - **Limited Functionality:** Callers expecting full ELSI functionality (like actual solutions to eigenvalue problems from `elsi_dm_real_sparse`) will not get meaningful results from those specific routines. However, the Fermi level and occupation calculation (`elsi_compute_mu_and_occ`) is functional and allows parts of the SCF cycle related to smearing to work as intended.
    - **`esl_smear_m`:** This module likely relies on `elsi_compute_mu_and_occ`, `elsi_set_mu_broaden_scheme`, and `elsi_set_mu_broaden_width` (or wrapper routines that call them) to determine electronic occupations and the Fermi energy.
    - **Build System:** The decision to use this fake module versus the real ELSI library is typically handled at the build system level (e.g., via preprocessor flags or by linking different object files).
