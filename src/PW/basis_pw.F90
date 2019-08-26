module esl_basis_pw_m
  use prec
  use yaml_output
  use fdf
  use esl_constants_m
  use esl_geometry_m
  use esl_utils_pw_m
  use esl_grid_m
  use esl_basis_base_m
 
  implicit none 
  private

  public :: basis_pw_t

  type, extends(basis_base_t) :: basis_pw_t
    real(dp) :: ecut !< Plane wave cut-off in Hartree
    integer  :: ndims(3) !< Number of plane-waves in each direction
    real(kind=dp) :: gmet(3,3) !< Metric

    real(kind=dp), allocatable :: gmod2(:) !length of the G-vectors
    integer, allocatable :: gmap(:,:) !Mapping
  contains
    private
    procedure, public :: init
    procedure, public :: summary
    final :: cleanup
  end type basis_pw_t
  
contains

  !Initialize the basis
  !----------------------------------------------------
  subroutine init(this, geo)
    class(basis_pw_t) :: this
    type(geometry_t), intent(in) :: geo

    integer :: i

    this%ecut = fdf_get('cut-off', 10._dp, 'Ha')

    ! Initialize auxiliary grid
    call ndims_from_ecut(this%ndims, this%ecut, geo%icell, [0._dp, 0._dp, 0._dp])
    call this%grid%init(this%ndims, geo%cell)

    ! Initialize PW basis
    do i = 1, 3
       this%gmet(i, :) = geo%icell(1, i) * geo%icell(1, :) &
              &        + geo%icell(2, i) * geo%icell(2, :) &
              &        + geo%icell(3, i) * geo%icell(3, :)
    end do

    this%size = get_number_of_pw(this%ndims, this%ecut, this%gmet, [0._dp, 0._dp, 0._dp])
 
    !TODO: We should create of these for each k-point
    allocate(this%gmod2(1:this%size))
    allocate(this%gmap(1:3,1:this%size))
    call construct_mod_map_tables(this%ndims, this%ecut, this%gmet, [0._dp, 0._dp, 0._dp], this%gmod2, this%gmap)

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(basis_pw_t) :: this

    if(allocated(this%gmod2)) deallocate(this%gmod2)
    if(allocated(this%gmap)) deallocate(this%gmap)

  end subroutine cleanup

  !Summary
  !----------------------------------------------------
  subroutine summary(this)
    class(basis_pw_t), intent(in) :: this

    call yaml_mapping_open("basis")
    call yaml_map("Type", "Plane waves")
    call yaml_map("Cut-off (Ha)", this%ecut)
    call yaml_map("Number of plane-waves", this%size)
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

end module esl_basis_pw_m
