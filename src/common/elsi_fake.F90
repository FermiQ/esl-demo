module elsi ! (Stub)
  use prec, only: ip, dp

  type elsi_handle
    integer :: tut
  end type elsi_handle

contains

  subroutine elsi_init(eh, solver, mode, mat_format, nbasis, nel, nstate)
    type(elsi_handle), intent(inout) :: eh
    real(dp), intent(in) :: nel
    integer(ip), intent(in) :: solver, mode, mat_format, nbasis, nstate
  end subroutine elsi_init

  subroutine elsi_finalize(eh)
    type(elsi_handle), intent(inout) :: eh
  end subroutine elsi_finalize

  subroutine elsi_set_csc(eh, global_nnz, local_nnz, nrow, col_idx, row_ptr)
    type(elsi_handle), intent(inout) :: eh
    integer(ip), intent(in) :: global_nnz, local_nnz, nrow
    integer(ip), intent(in) :: col_idx(:), row_ptr(:)
  end subroutine elsi_set_csc

  subroutine elsi_dm_real_sparse(eh, h, s, d, energy)
    type(elsi_handle), intent(in) :: eh
    real(dp), intent(inout) :: h(:), s(:), d(:)
    real(dp), intent(out) :: energy
  end subroutine elsi_dm_real_sparse

  subroutine elsi_get_edm_real_sparse(eh, edm)
    type(elsi_handle), intent(in) :: eh
    real(dp), intent(inout) :: edm(:)
  end subroutine elsi_get_edm_real_sparse

  subroutine elsi_compute_mu_and_occ(eh, nel, nstate, nspin, nkpt, evals, &
    & occ_nums, k_weights, fermi_level)
    type(elsi_handle), intent(in) :: eh
    real(dp), intent(in) :: nel
    integer(ip), intent(in) :: nstate, nspin, nkpt
    real(dp), dimension(nstate,nspin,nkpt), intent(in) :: evals
    real(dp), dimension(nstate,nspin,nkpt), intent(out) :: occ_nums
    real(dp), dimension(nkpt), intent(in) :: k_weights
    real(dp), intent(out) :: fermi_level
  end subroutine elsi_compute_mu_and_occ

  subroutine elsi_set_mu_broaden_scheme(eh, method)
    type(elsi_handle), intent(in) :: eh
    integer(ip), intent(in) :: method
  end subroutine elsi_set_mu_broaden_scheme

  subroutine elsi_set_mu_broaden_width(eh, width)
    type(elsi_handle), intent(in) :: eh
    real(dp), intent(in) :: width
  end subroutine elsi_set_mu_broaden_width

end module elsi
