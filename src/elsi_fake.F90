module elsi
  use prec
  type elsi_handle
     integer :: tut
  end type elsi_handle

  type sparse_matrix
     integer :: dat
  end type sparse_matrix

contains
  subroutine elsi_finalize(eh)
    type(elsi_handle), intent(in) :: eh
  end subroutine elsi_finalize
  subroutine elsi_dm_real(eh,d1, d2,d3,energy)
    type(elsi_handle), intent(in) :: eh
    integer, intent(inout) :: d1, d2, d3
    real(dp), intent(out) :: energy
  end subroutine elsi_dm_real
  subroutine elsi_compute_mu_and_occ(eh, n_electron, nstates, nspin, nkpt, &
       & eigenvalues, occ_numbers, k_weights, fermi_level)
    type(elsi_handle), intent(in) :: eh
     integer, intent(in) :: n_electron, nstates, nspin, nkpt
     real(dp), dimension(nstates,nspin,nkpt), intent(in) :: eigenvalues
     real(dp), dimension(nstates,nspin,nkpt), intent(in) :: occ_numbers
     real(dp), dimension(nkpt), intent(in) :: k_weights
     real(dp), intent(out) :: fermi_level
    
  end subroutine elsi_compute_mu_and_occ
end module elsi
