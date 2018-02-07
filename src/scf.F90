module scf_esl

  implicit none
  private

  public ::                       &
            scf_init,             &
            scf_end,              &
            scf_loop

  contains,
 

 !Initialize density and wfn
 !----------------------------------------------------
 subroutine scf_init( )


 end subroutine scf_init

 !Cleaning up
 !----------------------------------------------------
 subroutine scf_end( )

 end subroutine scf_end

 !Perform the self-consistent field calculation
 !----------------------------------------------------
 subroutine scf_loop( )

  
  integer :: max_iter  !< Maximum number of iterations
  integer :: iter !< Interation

  do iter = 1, max_iter
    !Diagonalization (ELSI/KSsolver)

    !Update occupations

    !Calc. density

    !Test tolerance and print status

    !Mixing (BLAS/LAPACK)

    !Update Hamiltonian matrix

  end do

 end subroutine scf_loop

end module scf_esl
