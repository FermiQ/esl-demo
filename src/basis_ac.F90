module esl_basis_ac_m
  use prec
  use yaml_output
  use esl_species_m

  implicit none

  private

  public :: basis_ac_t

  type basis_ac_t
    integer :: n_sites !< Number of different atomic centered sites
    integer :: n_species !< Number of different atomic centered species used on the sites
    integer :: n_functions !< Number of functions in basis (sum of number of functions per site)
    real(dp),        allocatable :: sites_xyz(:,:) !< Cartesian coordinates of the atomic centered sites
    type(species_t), allocatable :: species(:) !< The unique species used on the sites
    integer,         allocatable :: species_idx(:) !< A list of specie indices for each site
    integer,         allocatable :: site_function_start(:) !< Look-up table to convert a site index to the first global function index
    integer,         allocatable :: function_g2l(:) !< Look-up table to convert a global function index to the local function index on the specie
    integer,         allocatable :: function_site(:) !< Look-up table to convert a global function index to a site index.
  contains
    private
    procedure, public :: init
    procedure, public :: summary
    final :: cleanup
  end type basis_ac_t

contains

  !Initialize the basis
  !----------------------------------------------------
  subroutine init(this)
    class(basis_ac_t) :: this

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(basis_ac_t) :: this

  end subroutine cleanup

  !Summary
  !----------------------------------------------------
  subroutine summary(this)
    class(basis_ac_t), intent(in) :: this

    call yaml_mapping_open("basis_ac")
    call yaml_map("Number of functions", this%n_functions)
    call yaml_mapping_close()

  end subroutine summary

end module esl_basis_ac_m
