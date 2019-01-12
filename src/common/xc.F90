module esl_xc_m
  use prec, only : dp
  use esl_energy_m, only : energy_t
  use esl_geometry_m, only: geometry_t
  use esl_grid_m, only: grid_t

  use fdf, only : fdf_get

  ! GridXC initialization
  use gridxc, only: gridxc_init
  use gridxc, only: gridxc_setXC
  use gridxc, only: gridxc_cellXC
  use mesh3D, only: setMeshDistr  ! Defines a new mesh distribution
  use mesh3D, only: fftMeshDistr  ! Sets/gets distribution for FFTs
  use mesh3D, only: myMeshBox     ! Sets/gets distribution for FFTs
  
  use xc_f90_lib_m, only : XC_LDA_X, XC_LDA_C_PZ
  use xcmod, only : setxc_libxc_ids
#ifdef WITH_MPI
  use mpi, only: MPI_COMM_WORLD
#endif

  implicit none
  private

  public :: xc_t

  !Data structure for the xc potential
  type xc_t
    integer :: exchange
    integer :: correlation
    real(dp), pointer :: cell(:,:) => null()
    integer, pointer :: nmesh(:) => null()
  contains
    procedure, public :: init
    procedure, public :: calculate
    procedure, public :: summary
    final  :: cleanup
  end type xc_t

contains

  !Initialize the xc 
  !----------------------------------------------------
  subroutine init(this, geo, grid)
    class(xc_t), intent(inout) :: this
    type(geometry_t), target, intent(in) :: geo
    type(grid_t), target, intent(in) :: grid

    this%exchange = fdf_get('Exchange', XC_LDA_X)
    this%correlation = fdf_get('Correlation', XC_LDA_C_PZ)

    call setxc_libxc_ids(2, [this%exchange, this%correlation])
    
#ifdef WITH_MPI
    call gridxc_init(mpi_comm_world)
#else
    call gridxc_init()
#endif
    ! Specify flavour (FORCED TO LDA+CA)
    call gridxc_setXC(0, (/'LDA'/), (/'CA'/), (/1._dp/), (/1._dp/))

    this%cell => geo%cell
    this%nmesh => grid%ndims

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(xc_t), intent(inout) :: this

    this%cell => null()
    this%nmesh => null()

  end subroutine cleanup

  !Calc the xc potential from the density
  !----------------------------------------------------
  subroutine calculate(this, density, edata, vxc)
    class(xc_t), intent(in)  :: this
    real(dp),    intent(in)  :: density(:)
    real(dp),    intent(out) :: vxc(:)
    type(energy_t), intent(inout) :: edata

    integer  :: nspin
    integer  :: lb1, lb2, lb3, ub1, ub2, ub3
    integer  :: n_mesh(3)
    real(dp) :: ex, ec, dx, dc, stress_xc(3,3)

    ! FIXME: convert density from 1 to 4 dimensions

    ! Prepare bounding box explicitly (will be passed as argument when using MPI)
    ! We will read density(lb1:ub1, lb2:ub2, lb3:ub3, nspin)
    lb1 = 1
    ub1 = this%nmesh(1)
    lb2 = 1
    ub2 = this%nmesh(2)
    lb3 = 1
    ub3 = this%nmesh(3)

    ! Only valid for serial, will have to be passed when using MPI
    n_mesh(:) = [ub1, ub2, ub3]

    ! Set spin from last density dimension
    nspin = 1

    ! Let GridXC build the potential from the density
    call gridxc_cellXC( 1, this%cell, this%nmesh, &
        lb1, ub1, lb2, ub2, lb3, ub3, nspin, &
        density, Ex, Ec, Dx, Dc, stress_xc, Vxc)

    ! Populate energy terms, with necessary arithmetics
    ! Internally libgridxc uses Ry, convert to Hartree
    edata%exchange = ex * 0.5_dp
    edata%correlation = ec * 0.5_dp
    edata%int_nvxc = (dx + dc) * 0.5_dp
    vxc = vxc * 0.5_dp

    ! FIXME: convert potential from 4 to 1 dimension

  end subroutine calculate

  subroutine summary(this)
    use yaml_output
    class(xc_t), intent(in) :: this

    call yaml_mapping_open("XC")
    call yaml_map("Exchange (LibXC ID)", this%exchange)
    call yaml_map("Correlation (LibXC ID)", this%correlation)
    call yaml_mapping_close()

  end subroutine summary

end module esl_xc_m
