!> This module interfaces with ELSI to compute (currently) the density matrix,
!> and energy-weighted density matrix.
module esl_calc_density_matrix_ac_m

  use prec, only: ip, dp
  use esl_elsi_m, only: elsi_t
  use esl_sparse_matrix_m, only: sparse_matrix_t

#ifdef WITH_MPI
  use elsi, only: elsi_set_mpi
#endif

  use elsi, only: elsi_set_csc, &
                  elsi_dm_real_sparse, &
                  elsi_get_edm_real_sparse

  implicit none
  private

#ifdef WITH_MPI
  public :: set_elsi_sparsity_pattern_ac
#endif

  public :: calc_density_matrix_ac
  public :: calc_energy_density_matrix_ac

contains

#ifdef WITH_MPI
  subroutine set_elsi_sparsity_pattern_ac(elsi_ac, mpi_dist_ac, &
    & sparsity_pattern_ac)
    use mpi_dist, only: mpi_dist_t
    use esl_sparse_pattern_m, only: sparse_pattern_t

    type(elsi_t),           intent(inout) :: elsi_ac
    type(mpi_dist_t),       intent(in)    :: mpi_dist_ac
    type(sparse_pattern_t), intent(in)    :: sparsity_pattern_

    !Initialize ELSI MPI
    call elsi_set_mpi(elsi_ac%e_h, mpi_dist_ac%comm)

    !Initialize ELSI sparsity pattern
    call elsi_set_csc(elsi_ac%e_h, mpi_dist_ac%global_N, sparsity_pattern_%nz, &
      & sparsity_pattern_%nr, sparsity_pattern_%column, sparsity_pattern_%rptr)

  end subroutine set_elsi_sparsity_pattern_ac
#endif

  !Compute density matrix
  !----------------------------------------------------
  subroutine calc_density_matrix_ac(elsi_ac, mat_H, mat_S, mat_DM, KS_energy)

    type(elsi_t),          intent(inout) :: elsi_ac
    type(sparse_matrix_t), intent(inout) :: mat_H
    type(sparse_matrix_t), intent(inout) :: mat_S
    type(sparse_matrix_t), intent(inout) :: mat_DM
    real(dp),              intent(out)   :: KS_energy

    !Only real for now
    call elsi_dm_real_sparse(elsi_ac%e_h, mat_H%M, mat_S%M, mat_DM%M, KS_energy)

  end subroutine calc_density_matrix_ac

  !Compute energy-weighted density matrix
  !----------------------------------------------------
  subroutine calc_energy_density_matrix_ac(elsi_ac, mat_EDM)

    type(elsi_t),          intent(inout) :: elsi_ac
    type(sparse_matrix_t), intent(inout) :: mat_EDM

    !Only real for now
    call elsi_get_edm_real_sparse(elsi_ac%e_h, mat_EDM%M)

  end subroutine calc_energy_density_matrix_ac

end module esl_calc_density_matrix_ac_m
