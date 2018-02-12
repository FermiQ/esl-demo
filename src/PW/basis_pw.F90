module esl_basis_pw_m

  use prec
  use yaml_output
  use esl_geometry_m
  use esl_utils_pw_m
  use esl_grid_m
  
  private

  public :: basis_pw_t

  type basis_pw_t
    real(dp) :: ecut !< Plane wave cut-off in Hartree
    integer  :: npw  !< Number of plane waves
    integer  :: ndims(3) !< Number of plane-waves in each direction
    real(kind=dp) :: gmet(3,3) !< Metric
    integer*8 fftplan !< Forward FFT plan
    integer*8 ifftplan !< Backward FFT (IFFT) plan

    type(grid_t), pointer :: grid
  contains
    private
    procedure, public :: init
    procedure, public :: summary
    final :: cleanup
  end type basis_pw_t
  
contains

  !Initialize the basis
  !----------------------------------------------------
  subroutine init(this, grid, ecut, ndims, gcell)
    class(basis_pw_t)                :: this
    type(grid_t), target, intent(in) :: grid
    real(dp),             intent(in) :: ecut
    integer,              intent(in) :: ndims(3)
    real(dp),             intent(in) :: gcell(3,3)

    complex(dp),         allocatable :: arr(:)

    this%ecut = ecut

    do i = 1, 3
       this%gmet(i, :) = gcell(1, i) * gcell(1, :) &
              &        + gcell(2, i) * gcell(2, :) &
              &        + gcell(3, i) * gcell(3, :)
    end do


    this%npw = get_number_of_pw(ndims, ecut, this%gmet, [0._dp, 0._dp, 0._dp])

    this%ndims(1:3) = ndims(1:3)

    this%grid => grid

    ! Initialization for FFT and IFFT
    ! TODO include 'fftw3.f90' should be put properly
    allocate(arr(ndims(1),ndims(2),ndims(3)))

    call dfftw_plan_dft_3d(this%fftplan, ndims(1), ndims(2), ndims(3), &
      arr, arr, FFTW_FORWARD, FFTW_ESTIMATE)

    call dfftw_plan_dft_3d(this%ifftplan, ndims(1), ndims(2), ndims(3), &
      arr, arr, FFTW_BACKWARD, FFTW_ESTIMATE)

    deallocate(arr)

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(basis_pw_t) :: this
    ! Deconstructor for fft plan
    call dfftw_destroy_plan(this%fftplan)
    call dfftw_destroy_plan(this%ifftplan)
  end subroutine cleanup

  !Summary
  !----------------------------------------------------
  subroutine summary(this)
    class(basis_pw_t), intent(in) :: this

    call yaml_mapping_open("basis_pw")
    call yaml_map("Cut-off (Ha)", this%ecut)
    call yaml_map("Number of plane-waves", this%npw)
    call yaml_mapping_close()

  end subroutine summary

end module esl_basis_pw_m
