module basis_esl
  use fdf, only : block_fdf, fdf_block,fdf_defined, &
                  parsed_line, fdf_breals, fdf_bline, fdf_bnames, &
                  fdf_get

  use message_esl

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

     character(len=100) :: str
     str=fdf_get('BasisSet','Planewaves')
     if (str == 'Planewaves')  Then
       this%basis_type=PLANEWAVES
     else if (str=='AtomicOrbitals') then
       this%basis_type=ATOMICORBS
     else
       call message_error("The specified basis set is not correct.")
     endif

   end subroutine init
 
   !Release
   !----------------------------------------------------
   subroutine cleanup(this)
     type(basis_t) :: this

   end subroutine cleanup

end module basis_esl
