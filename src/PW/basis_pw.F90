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
    real(kind=dp), allocatable :: gmod(:) !length of the G-vectors
    integer, allocatable :: gmap(:,:) !Mapping
 
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

    this%ecut = ecut

    do i = 1, 3
       this%gmet(i, :) = gcell(1, i) * gcell(1, :) &
              &        + gcell(2, i) * gcell(2, :) &
              &        + gcell(3, i) * gcell(3, :)
    end do

    this%npw = get_number_of_pw(ndims, ecut, this%gmet, [0._dp, 0._dp, 0._dp])
    
    !TODO: We should create of these for each k-point
    allocate(this%gmod(1:npw))
    allocate(this%gmap(1:3,1:npw))
    call construct_mod_map_tables(ndims, ecut, this%gmet, [0._dp, 0._dp, 0._dp], this%gmod, this%gmap)

    this%ndims(1:3) = ndims(1:3)

    this%grid => grid

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(basis_pw_t) :: this

    if(allocated(this%gmod)) deallocate(this%gmod)
    if(allocated(this%gmap)) deallocate(this%gmap)

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
