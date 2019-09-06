module esl_basis_ac_abc_t
  use prec
  use pspiof_m

  use esl_basis_base_m, only: basis_base_t
  use esl_sparse_pattern_m, only: sparse_pattern_t
  use esl_sparse_matrix_m, only: sparse_matrix_t

  implicit none

  private

  public :: orbital_ac_t
  public :: state_ac_t
  public :: basis_ac_abc_t


  ! Local container for the orbital basis function.
  ! This orbital object contains three things:
  !  1. l, l-quantum number
  !  2. m, m-quantum number
  !  3. R(r), radial function.
  !  4. occ, initial occupation
  ! From these 3 quantities one can re-create the full
  !   psi(r) = R(r) Y_l^m(r)
  ! values at any point.
  ! Later this orbital may contain information such
  ! as polarization, zeta-information etc.
  type orbital_ac_t
    
    integer :: l = 0, m = 0
    real(dp) :: r_cut = 0._dp
    real(dp) :: occ = 0._dp
    !< Radial function of psi [R(r) Y_l^m(r) == psi(r)]
    type(pspiof_meshfunc_t), pointer :: R => null()
    
  end type orbital_ac_t
  
  ! Each site consists of a set of orbitals
  ! These orbitals are the *expanded* versions of a state
  ! I.e. n_orbital is the actual total number of basis-orbitals on that
  ! site.
  type state_ac_t

    ! Number of orbitals
    integer :: n_orbital = 0
    ! Maximum r_cut for all orbitals on this state
    real(dp) :: r_cut = 0._dp
    type(orbital_ac_t), pointer :: orb(:) => null()
    
  end type state_ac_t

  ! In order to call methods with circular dependencies we *have* to
  ! have an abstract class definition outside of the actual usage
  ! Probably this should be carefully thought out to really have a good
  ! class dependency structure.
  ! This should probably contain procedure names such that we force
  ! certain procedures
  type, abstract, extends(basis_base_t) :: basis_ac_abc_t

    ! TODO this should be from system, i.e.
    real(dp) :: Q = 0._dp !< Number of electrons

    integer :: n_site = 0 !< Number of atomic centered sites
    real(dp), allocatable :: xyz(:,:) !< Cartesian coordinates of the atomic centered sites

    !< Number of unique sites (irrespective of Cartesian positions)
    integer :: n_state = 0
    !< Container for each site. Each site contains a set of orbitals on which the basis is expanded
    type(state_ac_t), allocatable :: state(:)

    !< Total number of orbitals in this AC-basis
    !< This equates to sum[ state(is)%n_orbital * <number of sites with state is> for all is]
    integer :: n_orbital = 0 !< Number of functions in basis (sum of number of functions per site)

    integer, allocatable :: site_state_idx(:) !< State indices for each site

    integer, allocatable :: site_orbital_start(:) !< Look-up table to convert a site index to the first global orbital index (n_site + 1)
    integer, allocatable :: orbital_site(:) !< Look-up table to convert a global orbital index to a site index.

    ! Currently we don't need this array.
    ! It may be used later.
    !    integer, allocatable :: orbital_g2l(:) !< Look-up table to convert a global orbital index to the local function index on the specie

    !< Sparsity pattern describing where we have non-zero elements between orbitals
    type(sparse_pattern_t) :: sparse_pattern
    !< Non-orthogonal basis (overlap matrix) is defined uniquely from the basis (and geometry!)
    type(sparse_matrix_t) :: S ! always 1D

  end type basis_ac_abc_t
  
end module esl_basis_ac_abc_t
