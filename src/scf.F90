module scf_esl
  use prec, only : dp,ip

  implicit none
  private

  public ::                       &
            scf_t,                &
            scf_loop

  !Data structure containing the data for the SCF
  type scf_t 
   real(kind=dp) :: tol_reldens
   contains
     private
     procedure, public :: init
     final :: cleanup
  end type scf_t

  contains
 

 !Initialize density and wfn
 !----------------------------------------------------
 subroutine init(this)
   class(scf_t)  :: this

   !Parse here the data for the SCF
 end subroutine init

 !Cleaning up
 !----------------------------------------------------
 subroutine cleanup(this)
   type(scf_t) :: this


 end subroutine cleanup

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
