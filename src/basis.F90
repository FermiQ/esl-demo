module esl_basis_m

  use prec
  use fdf, only : fdf_get, leqi
  use message_esl
  use esl_basis_pw_m

  implicit none

  private

  public :: basis_t

  !Data structure for the basis
  type basis_t
     integer :: type
     type(basis_pw_t) :: pw
   contains
     private
     procedure, public :: init
     procedure, public :: summary
     final :: cleanup
  end type basis_t

  integer, public, parameter :: &
       PLANEWAVES   = 1, &
       ATOMCENTERED = 2

contains

  !Initialize the basis set
  !----------------------------------------------------
  subroutine init(this)
    class(basis_t) :: this

    character(len=100) :: str
    real(dp) :: ecut
    integer  :: ndims(3)
    real(dp) :: gcell(3,3)

    str = fdf_get('BasisSet', 'Planewaves')
    if ( leqi(str, 'Planewaves') ) then
       this%type = PLANEWAVES

       ecut = fdf_get('cut-off', 10._dp, 'Ha')
       call this%pw%init(ecut, ndims, gcell)

    else if ( leqi(str, 'AtomicOrbitals') ) then
       this%type = ATOMCENTERED

    else
       call message_error("Unknown basis set: "//trim(str))
    end if

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(basis_t) :: this

  end subroutine cleanup

  !Summary
  !----------------------------------------------------
  subroutine summary(basis)
    use yaml_output
    class(basis_t), intent(in) :: basis

    call yaml_mapping_open("basis")
    select case (basis%type)
    case ( PLANEWAVES )
       call yaml_map("Type", "Plane waves")
       call basis%pw%summary()
    case ( ATOMCENTERED )
       call yaml_map("Type", "Atomic orbitals")
    end select
    call yaml_mapping_close()

  end subroutine summary

end module esl_basis_m
