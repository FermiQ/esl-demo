module scf_esl
  use parser_esl

  implicit none
  private

  public ::                       &
            scf_t,                &
            scf_init,             &
            scf_end,              &
            scf_loop

  !Data structure containing the data for the SCF
  type scf_t 
   double precision :: tol_reldens
  end type scf_t

  contains
 

 !Initialize density and wfn
 !----------------------------------------------------
 subroutine scf_init(this, parser)
   type(scf_t),  intent(inout) :: this
   type(parser_t),  intent(in) :: parser

   !Parse here the data for the SCF
 end subroutine scf_init

 !Cleaning up
 !----------------------------------------------------
 subroutine scf_end(this)
   type(scf_t),  intent(inout) :: this


 end subroutine scf_end

 !Perform the self-consistent field calculation
 !----------------------------------------------------
 subroutine scf_loop(this)
   type(scf_t),  intent(inout) :: this

  
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
