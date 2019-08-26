module esl_basis_m
  use prec
  use fdf, only : fdf_get, leqi
  use esl_message_m
  use esl_basis_ac_m
  use esl_basis_pw_m
  use esl_constants_m
  use esl_geometry_m
  use esl_grid_m
  use yaml_output

#ifdef WITH_FLOOK
  use dictionary
  use esl_dict_m
#endif

  implicit none

  private

  public :: basis_t

  !< Basis information
  !<
  !< Contains the choice of the used basis, PW vs. AC
  !< and also contains the specific basises used for each
  !< of the PW vs. AC codes.
  type basis_t
    integer :: type
    type(basis_pw_t) :: pw !< Plane-wave basis
    type(basis_ac_t) :: ac !< Atomic-centered basis
    type(grid_t)     :: grid !< Auxiliary real-space grid
    integer :: size
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
  subroutine init(this, geo)
    class(basis_t)                   :: this
    type(geometry_t),     intent(in) :: geo
    
    character(len=128) :: str
    real(dp) :: ecut
    integer  :: ndims(3)

    str = fdf_get('BasisSet', 'Planewaves')
    if ( leqi(str, 'Planewaves') ) then
      
      this%type = PLANEWAVES

      ! Initialize PW basis
      call this%pw%init(geo)

#ifdef WITH_FLOOK
      call esl_dict_var_add('PW.CutOff', this%pw%ecut)
#endif
      
    else if ( leqi(str, 'AtomicOrbitals') ) then
      
      this%type = ATOMCENTERED

      ! Initialize AC basis
      call this%ac%init(geo)
      
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
  subroutine summary(this)
    class(basis_t), intent(in) :: this

    select case (this%type)
    case ( PLANEWAVES )
      call this%pw%summary()
    case ( ATOMCENTERED )
      call this%ac%summary()
    end select

  end subroutine summary

  subroutine ndims_from_spacing(ndims, spacing, cell)
    integer,  intent(out) :: ndims(3)
    real(dp), intent(in)  :: spacing
    real(dp), intent(in)  :: cell(3,3)

    integer :: idim, n
    
    do idim = 1,3
      n = ceiling(cell(idim, idim) / spacing) - 1
      call fourier_dim(n, ndims(idim))
    end do

  end subroutine ndims_from_spacing
    
end module esl_basis_m
