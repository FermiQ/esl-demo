module basis_esl

  implicit none
  private

  public ::                &
            basis_t
  
  !Data structure for the basis
  type basis_t
    integer :: basis_type
    integer :: size !< Number of coefficients in the basis
    contains
      private
      procedure, public :: init
      final :: cleanup
  end type basis_t

  integer, public, parameter :: &
    PLANEWAVES   = 1,           &
    ATOMICORBS   = 2


  contains

   !Initialize the physical system
   !----------------------------------------------------
   subroutine init(this)
     class(basis_t) :: this

     !Parse the informations from the input file

   end subroutine init
 
   !Release
   !----------------------------------------------------
   subroutine cleanup(this)
     type(basis_t) :: this

   end subroutine cleanup

end module basis_esl
