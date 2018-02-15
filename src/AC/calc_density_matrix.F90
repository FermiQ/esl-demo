! This module interfaces with ELSI to compute (currently) the density matrix,
!> and energy-weighted density matrix.
!> TODO: Unify naming scheme.
module esl_calc_density_matrix_ac_m

  use prec, only: ip, dp
  use esl_elsi_m, only: elsi_t
  use esl_sparse_matrix_m, only: sparse_matrix_t

#ifdef WITH_MPI
  use elsi, only: elsi_set_mpi
#endif

  use elsi, only: elsi_set_output, &
                  elsi_set_csc, &
                  elsi_set_csc_blk, &
                  elsi_dm_real_sparse, &
                  elsi_get_edm_real_sparse, &
                  elsi_get_mu, &
                  elsi_get_entropy

  implicit none
  private

#ifdef WITH_MPI
  public :: set_elsi_sparsity_pattern_ac
#endif

  public :: calc_density_matrix_ac
  public :: calc_energy_density_matrix_ac
  public :: get_energy_results

contains

#ifdef WITH_MPI
  subroutine set_elsi_sparsity_pattern_ac(elsi_ac, mpi_dist_bc_ac, &
    & sparsity_pattern_ac)
    use mpi_dist_block_cyclic_m, only: mpi_dist_block_cyclic_t
    use esl_sparse_pattern_m, only: sparse_pattern_t

    type(elsi_t),                  intent(inout) :: elsi_ac
    type(mpi_dist_block_cyclic_t), intent(in)    :: mpi_dist_bc_ac
    type(sparse_pattern_t),        intent(in)    :: sparsity_pattern_ac

    ! Initialize ELSI MPI
    call elsi_set_mpi(elsi_ac%e_h, mpi_dist_bc_ac%comm)

    ! Initialize ELSI sparsity pattern
    call elsi_set_csc(elsi_ac%e_h, mpi_dist_bc_ac%global_N, &
      & sparsity_pattern_ac%nz, sparsity_pattern_ac%nr, &
      & sparsity_pattern_ac%column, sparsity_pattern_ac%rptr)

    ! Set block size
    call elsi_set_csc_blk(elsi_ac%e_h, mpi_dist_bc_ac%block)

    ! Request debug info
    if (mpi_dist_bc_ac%rank == 0) then
      call elsi_set_output(elsi_ac%e_h, 3)
    end if

  end subroutine set_elsi_sparsity_pattern_ac
#endif

  !> Compute density matrix
  subroutine calc_density_matrix_ac(elsi_ac, mat_H, mat_S, mat_DM)

    type(elsi_t),          intent(inout) :: elsi_ac
    type(sparse_matrix_t), intent(inout) :: mat_H
    type(sparse_matrix_t), intent(inout) :: mat_S
    type(sparse_matrix_t), intent(inout) :: mat_DM

    !Only real for now
    call elsi_dm_real_sparse(elsi_ac%e_h, mat_H%M, mat_S%M, mat_DM%M, &
      & elsi_ac%KS_energy)

  end subroutine calc_density_matrix_ac

  !> Compute energy-weighted density matrix
  subroutine calc_energy_density_matrix_ac(elsi_ac, mat_EDM)

    type(elsi_t),          intent(inout) :: elsi_ac
    type(sparse_matrix_t), intent(inout) :: mat_EDM

    !Only real for now
    call elsi_get_edm_real_sparse(elsi_ac%e_h, mat_EDM%M)

  end subroutine calc_energy_density_matrix_ac

  !> Retrive energies
  subroutine get_energy_results(elsi_ac, KS_energy, fermi_level, entropy)

    type(elsi_t), intent(inout) :: elsi_ac
    real(dp),     intent(out)   :: KS_energy
    real(dp),     intent(out)   :: fermi_level
    real(dp),     intent(out)   :: entropy

    call elsi_get_mu(elsi_ac%e_h, elsi_ac%fermi_level)
    call elsi_get_entropy(elsi_ac%e_h, elsi_ac%entropy)

    KS_energy   = elsi_ac%KS_energy
    fermi_level = elsi_ac%fermi_level
    entropy     = elsi_ac%entropy

  end subroutine get_energy_results

end module esl_calc_density_matrix_ac_m
