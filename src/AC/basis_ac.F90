module esl_basis_ac_m
  use prec
  use yaml_output

  use esl_message_m
  use esl_species_m
  use esl_geometry_m
  use esl_sparse_matrix_m
  
  implicit none

  private

  public :: basis_ac_t

  type basis_ac_t
    integer :: n_sites !< Number of different atomic centered sites
    integer :: n_species !< Number of different atomic centered species used on the sites
    integer :: n_functions !< Number of functions in basis (sum of number of functions per site)
    real(dp),        allocatable :: xyz(:,:) !< Cartesian coordinates of the atomic centered sites
    type(species_t), allocatable :: species(:) !< The unique species used on the sites
    integer,         allocatable :: species_idx(:) !< A list of specie indices for each site
    integer,         allocatable :: site_function_start(:) !< Look-up table to convert a site index to the first global function index (n_sites + 1)
    integer,         allocatable :: function_g2l(:) !< Look-up table to convert a global function index to the local function index on the specie
    integer,         allocatable :: function_site(:) !< Look-up table to convert a global function index to a site index.
  contains
    private
    procedure, public :: init
    procedure, public :: summary

    procedure, public :: atomic_density_matrix

    final :: cleanup
  end type basis_ac_t

contains

  !Initialize the basis
  !----------------------------------------------------
  subroutine init(this, geo)
    class(basis_ac_t) :: this
    type(geometry_t), intent(in) :: geo

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(basis_ac_t) :: this

  end subroutine cleanup

  !< Create a sparse DM containing the atomic fillings
  !<
  !< This atomic filling density matrix can easily be used
  !< to calculate the neutral atom density as well as using it
  !< to subtract the neutral atom density from the full
  !< density matrix.
  !< I.e.
  !<
  !< \code{.f90}
  !<  call density_ac%add_density_matrix(grid, basis, DM)
  !<  call basis%atomic_density_matrix(atomic_DM)
  !<  atomic_DM%M = - atomic_DM%M
  !<  call density_ac%add_density_matrix(grid, basis, atomic_DM)
  !<  ! at this point rho contains the dRho = Rho - AtomicRho
  !< \endcode
  !<
  !< NOTE!
  !< If DM has not been initialized on entry, the sparse pattern that
  !< gets associated with this DM is the
  !< only reference in the program. I.e. when this DM is not to
  !< be used anymore the DM%sp object has to be deleted manually.
  subroutine atomic_density_matrix(this, DM)
    class(basis_ac_t), intent(in) :: this
    type(sparse_matrix_t), intent(inout) :: DM

    ! Local variables
    type(sparse_pattern_t), pointer :: sp => null()

    integer :: io, ind

    if ( .not. DM%initialized() ) then
      
      ! This routine call will create an initial DM matrix with
      ! minimal sparse pattern. It should only be used
      ! for cases when the atomic DM is the only thing to worry
      ! about. I.e. subsequent usage of the sparse matrix is
      ! ill-adviced for anything but initial Mulliken, initial
      ! density on grids etc.
      call init()
      
    else
      
      call message_error('atomic_density_matrix: non-initialized DM &
          &is Currently not implemented!')
      
    end if

    ! Initialize everything to 0
    DM%M(:) = 0._dp
    
    ! Fill the density matrix with 1's
    do io = 1, this%n_functions
      
      ! Figure out if this basis function has initial density
      do ind = sp%rptr(io) , sp%rptr(io) + sp%nrow(io) - 1

        ! Check for the diagonal part
        if ( sp%column(ind) == io ) then
          
          ! Ok, we have a diagonal entry.
          DM%M(ind) = 1._dp
          
        end if
        
      end do
      
    end do

    nullify(sp)
    
  contains
      
    subroutine init()
      integer :: ia

      ! Create the sparse pattern
      allocate(sp)

      ! Initialize the sparse pattern with only one element per-row
      call sp%init(this%n_functions, this%n_functions, np=1)

      ! Add elements to the sparse matrix
      do ia = 1, this%n_sites

        ! TODO fill here to correctly fill the diagonal elements
        ! we should loop on the electronic configuration

      end do

      ! Finalize the sparse object to remove all unnecessary elements
      call sp%finalize()

      ! Note that now the DM object contains the *only* reference
      ! to sp that exists.
      call DM%init(sp)

    end subroutine init
    
  end subroutine atomic_density_matrix


  !Summary
  !----------------------------------------------------
  subroutine summary(this)
    class(basis_ac_t), intent(in) :: this

    integer :: is
    character(len=10) :: str

    call yaml_mapping_open("basis_ac")
    
    call yaml_map("Number of sites", this%n_sites)
    
    call yaml_sequence_open("Site coordinates", advance = "no")
    call yaml_comment("| X | Y | Z |", hfill = "-")
    do is = 1, this%n_sites
      call yaml_sequence(advance="no")
      write(str, '(i0)') is
      call yaml_map(trim(str), this%xyz(:,is))
    end do
    call yaml_sequence_close()
    
    call yaml_map("Number of functions", this%n_functions)
    call yaml_mapping_close()

  end subroutine summary

end module esl_basis_ac_m
