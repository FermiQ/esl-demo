module basis_esl

  use prec
  use fdf, only : fdf_get, leqi

  use message_esl

  implicit none
  private

  public ::                &
       basis_t

  type pw_basis_t
     real(dp) :: ecut
  end type pw_basis_t

  !Data structure for the basis
  type basis_t
     integer :: type
     type(pw_basis_t) :: pw_basis
    integer :: size !< Number of coefficients in the basis
    contains
      private
      procedure, public :: init
      procedure, public :: init_atomic_orbitals
      procedure, public :: summary
      final :: cleanup
  end type basis_t

  integer, public, parameter :: &
       PLANEWAVES   = 1, &
       ATOMICORBS   = 2


  contains

   !Initialize the physical system
   !----------------------------------------------------
   subroutine init(this)
     class(basis_t) :: this

     character(len=100) :: str
     real(dp) :: ecut

     str = fdf_get('BasisSet','Planewaves')
     if ( leqi(str, 'Planewaves') )  Then
        this%type = PLANEWAVES
        this%pw_basis%ecut = fdf_get('cut-off', 10._dp, 'Ha')
     else if ( leqi(str, 'AtomicOrbitals') ) then
        this%type = ATOMICORBS
     else
        call message_error("The specified basis set is not correct.")
     endif

     !TEMP
     this%size = 1

   end subroutine init
 
   !Release
   !----------------------------------------------------
   subroutine cleanup(this)
     type(basis_t) :: this

   end subroutine cleanup

   !Initialize the atomic orbitals
   !----------------------------------------------------
   subroutine init_atomic_orbitals(this)
     class(basis_t) :: this

     if( this%type /= ATOMICORBS ) return

   end subroutine init_atomic_orbitals
   
   !Summary
   !----------------------------------------------------
   subroutine summary(basis)
     use yaml_output
     class(basis_t), intent(in) :: basis

     call yaml_mapping_open("basis")
     select case (basis%type)
     case ( PLANEWAVES )
        call yaml_map("Type", "Plane waves")
        call yaml_map("Cut-off (Ha)", basis%pw_basis%ecut)
     case ( ATOMICORBS )
        call yaml_map("Type", "Atomic orbitals")
     end select
     call yaml_mapping_close()
   end subroutine summary

end module basis_esl
