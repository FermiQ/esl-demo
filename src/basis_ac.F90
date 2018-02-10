module esl_basis_ac_m
  use prec
  use yaml_output
  use esl_species_m

  implicit none
    
  private

  public :: basis_ac_t

  type basis_ac_t
    integer  :: n_sites
    integer  :: n_species
    integer  :: n_functions  !< Number of functions in basis
    real(dp),        allocatable :: sites_xyz(:,:)
    type(species_t), allocatable :: species(:)
    integer,         allocatable :: species_idx(:)
    integer,         allocatable :: site_function_start(:)
    integer,         allocatable :: function_g2l(:)
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
