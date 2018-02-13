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

  implicit none

  private

  public :: basis_t

  !Data structure for the basis
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
    
    character(len=100) :: str
    real(dp) :: ecut
    integer  :: ndims(3)

    str = fdf_get('BasisSet', 'Planewaves')
    if ( leqi(str, 'Planewaves') ) then
      
      this%type = PLANEWAVES
      ecut = fdf_get('cut-off', 10._dp, 'Ha')

      ! Initialize auxiliary grid
      call ndims_from_ecut(ndims, ecut, geo%icell, [0._dp, 0._dp, 0._dp])
      call this%grid%init(ndims, geo%cell)

      ! Initialize PW basis
      call this%pw%init(this%grid, ecut, ndims, geo%icell)
      this%size = this%pw%npw
      
    else if ( leqi(str, 'AtomicOrbitals') ) then
      
      this%type = ATOMCENTERED

      ! Initialize auxiliary grid
      ! For the moment the spacing in real space is hardcoded
      call ndims_from_spacing(ndims, 0.25_dp, geo%cell)
      call this%grid%init(ndims, geo%cell)

      ! Initialize AC basis
      call this%ac%init(geo)
      this%size = this%ac%n_orbital
      
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

    call yaml_mapping_open("basis")
    select case (this%type)
    case ( PLANEWAVES )
      call yaml_map("Type", "Plane waves")
      call this%pw%summary()
    case ( ATOMCENTERED )
      call yaml_map("Type", "Atomic orbitals")
      call this%ac%summary()
    end select
    call yaml_mapping_close()

    call this%grid%summary()
    
  end subroutine summary

  subroutine ndims_from_ecut(ndims, ecut, gcell, kpt)
    integer,  intent(out) :: ndims(3)
    real(dp), intent(in)  :: ecut
    real(dp), intent(in)  :: gcell(3,3)
    real(dp), intent(in)  :: kpt(3)

    real(dp) :: threshold
    integer  :: dir, i
    real(dp) :: gmet(3,3)
    real(dp), parameter :: BOXCUTMIN = 2.0_dp

    do i = 1, 3
      gmet(i, :) = gcell(1, i) * gcell(1, :) + &
                   gcell(2, i) * gcell(2, :) + &
                   gcell(3, i) * gcell(3, :)
    end do

    threshold = 0.5_dp * BOXCUTMIN**2 * ecut / PI**2
    !@todo: Don't take into account symmetries or k points.
    ndims = 16
    do
      if (smallest(ndims, gmet, dir) >= threshold) exit
      call fourier_dim(ndims(dir) + 1, ndims(dir))
    end do

  contains

    function smallest(ndims, gmet, dir)
      integer,  intent(in) :: ndims(3)
      real(dp), intent(in) :: gmet(3,3)
      integer,  intent(out) :: dir
      real(dp) :: smallest

      integer :: i1, i2, i3, s1, s2, s3, idir
      real(dp) :: prev

      smallest = dsq(ndims(1) / 2, -ndims(2) / 2, -ndims(3) / 2) + 0.01_dp
      do idir = 1, 3
        s1 = 0
        if (idir == 1) s1 = ndims(1)/2
        s2 = 0
        if (idir == 2) s2 = ndims(2)/2
        s3 = 0
        if (idir == 3) s3 = ndims(3)/2
        do i1 = s1, ndims(1)/2
          do i2 = s2, ndims(2)/2
            do i3 = s3, ndims(3)/2
              prev = smallest
              smallest = min(smallest, dsq( i1,  i2,  i3))
              smallest = min(smallest, dsq(-i1,  i2,  i3))
              smallest = min(smallest, dsq( i1, -i2,  i3))
              smallest = min(smallest, dsq(-i1, -i2,  i3))
              smallest = min(smallest, dsq( i1,  i2, -i3))
              smallest = min(smallest, dsq(-i1,  i2, -i3))
              smallest = min(smallest, dsq( i1, -i2, -i3))
              smallest = min(smallest, dsq(-i1, -i2, -i3))
              if (prev /= smallest) dir = idir
            end do
          end do
        end do
      end do
    end function smallest

    function dsq(i1, i2, i3)
      integer, intent(in) :: i1, i2, i3
      real(dp) :: dsq

      dsq = gmet(1,1)*(kpt(1) + dble(i1))**2 &
          + gmet(2,2)*(kpt(2) + dble(i2))**2 &
          + gmet(3,3)*(kpt(3) + dble(i3))**2 &
          + 2._dp*(gmet(1,2)*(kpt(1) + dble(i1))*(kpt(2) + dble(i2)) &
          + gmet(2,3)*(kpt(2) + dble(i2))*(kpt(3) + dble(i3)) &
          + gmet(3,1)*(kpt(3) + dble(i3))*(kpt(1) + dble(i1)))

    end function dsq
  end subroutine ndims_from_ecut
  
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
