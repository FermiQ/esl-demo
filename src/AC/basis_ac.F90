module esl_basis_ac_m
  use prec
  use yaml_output
  use pspiof_m

  use esl_numeric_m
  use esl_message_m
  use esl_species_m
  use esl_geometry_m
  use esl_sparse_matrix_m
  
  implicit none

  private

  public :: basis_ac_t


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
    real(dp) :: occ = 0._dp
    !< Radial function of psi [R(r) Y_l^m(r) == psi(r)]
    type(pspiof_meshfunc_t), pointer :: R => null()
    
  end type orbital_ac_t
  
  ! Each site consists of a set of orbitals
  ! These orbitals are the *expanded* versions of a state
  ! I.e. n_orbital is the actual total number of basis-orbitals on that
  ! site.
  type state_ac_t
    
    integer :: n_orbital = 0
    type(orbital_ac_t), pointer :: orb(:) => null()
    
  end type state_ac_t

  type basis_ac_t

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

  contains
    private
    procedure, public :: init
    procedure, public :: summary

    procedure, private :: get_psi_all, get_psi_single
    generic, public :: get_psi => get_psi_all, get_psi_single
    
    procedure, public :: atomic_density_matrix

    final :: cleanup
    
  end type basis_ac_t


contains

  !Initialize the basis
  !----------------------------------------------------
  subroutine init(this, geo)
    class(basis_ac_t) :: this
    type(geometry_t), intent(in) :: geo

    ! Local variables
    type(pspiof_meshfunc_t), pointer :: mesh_R => null()
    integer :: is, no, io, isite
    integer :: l, m
    real(dp) :: occ

    ! Construct the basis functions

    ! First copy over the Cartesian coordinates
    this%n_site = geo%n_atoms
    allocate( this%xyz(3, this%n_site) )
    
    ! Copy coordinates
    this%xyz(:,:) = geo%xyz(:,:)

    ! Site to state pointers
    ! This is equivalent to the species_idx in geo%
    allocate( this%site_state_idx(this%n_site) )
    this%site_state_idx = geo%species_idx

    ! Starting orbital pointers according to each site
    allocate( this%site_orbital_start(this%n_site + 1) )

    ! Each specie corresponds to a "state"
    this%n_state = geo%n_species
    allocate( this%state(this%n_state) )
    
    ! Loop over each specie and construct the basis orbitals
    do is = 1, geo%n_species
      
      ! For each state we need to figure out how many basis-orbitals
      ! the state is expanded into.
      ! So the species%state does not per see contain all orbitals

      no = 0
      ! %n_orbitals is a bit misleading as they are not orbitals.
      ! Rather, they are l, R states
      do io = 1, geo%species(is)%n_radial_orbitals

        call geo%species(is)%get_radial_orbital(io, ll=l)
        no = no + l * 2 + 1

      end do

      ! Allocate all orbitals for this site
      ! Note that this is aranged according to the species
      ! states
      !   for species%rad_orb has
      !     rad_orb , l = 1
      !     rad_orb , l = 2
      !     rad_orb , l = 1
      ! the final orbitals will be:
      !   l=1,m=-1
      !   l=1,m=0
      !   l=1,m=1
      !   l=2,m=-2
      !   l=2,m=-1
      !   l=2,m=0
      !   l=2,m=1
      !   l=2,m=2
      !   l=1,m=-1
      !   l=1,m=0
      !   l=1,m=1
      this%state(is)%n_orbital = no
      allocate( this%state(is)%orb(no) )

      ! Loop and create the things
      no = 0 ! counter for the current placement of the radial orbital
      do io = 1, geo%species(is)%n_radial_orbitals

        ! Retrieve it once, then subsequent orbitals are pointing to a single
        ! instance.
        allocate(mesh_r)

        ! Get all radial orbital information
        call geo%species(is)%get_radial_orbital(io, l, mesh_r, occ)

        do m = -l, l
          ! Populate the orbital
          no = no + 1
          this%state(is)%orb(no)%l = l
          this%state(is)%orb(no)%m = m
          this%state(is)%orb(no)%R => mesh_r
          this%state(is)%orb(no)%occ = occ / (l * 2 + 1)
        end do

        ! Nullify to retain only the references within the orbital lists
        nullify(mesh_r)

      end do
      
    end do

    ! Calculate total number of orbitals
    no = 0
    do isite = 1, this%n_site
      
      is = this%site_state_idx(isite)
      no = no + this%state(is)%n_orbital
      
    end do
          
    ! Allocate an orbital -> site index
    allocate( this%orbital_site(no) )
    
    ! Populate orbital indices
    no = 1
    do isite = 1, this%n_site

      this%site_orbital_start(isite) = no

      ! Accummulate
      is = this%site_state_idx(isite)
      do io = 1, this%state(is)%n_orbital
        this%orbital_site(no) = isite
        no = no + 1
      end do

    end do
    this%site_orbital_start(this%n_site+1) = no
    this%n_orbital = no - 1

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(basis_ac_t) :: this
    integer :: is, io
    type(pspiof_meshfunc_t), pointer :: R => null()

    if ( .not. allocated(this%xyz) ) return

    deallocate(this%xyz)
    deallocate(this%site_state_idx)
    deallocate(this%site_orbital_start)
    deallocate(this%orbital_site)

    ! Loop over all states
    do is = 1, this%n_state

      ! Loop all orbitals on this state
      do io = 1, this%state(is)%n_orbital

        ! This should find the unique orbitals and deallocate each of them
        ! Also see init
        ! In init we associate the same radial function %orb(..)%R with
        ! multiple orbitals, since for each l we have -l:l different m quantum
        ! numbers with the same radial function.
        ! Instead of duplicating each of them, we use pointers to
        ! create a single instance and let all orbitals point to the same
        ! thing. This complicates things a bit.
        ! One could also make an l-state-container, but that would increase
        ! the type-levels by 1.
        if ( .not. associated(R) ) then
          R => this%state(is)%orb(io)%R

          cycle
          
        else if ( associated(this%state(is)%orb(io)%R, R) ) then
          nullify( this%state(is)%orb(io)%R )
          
          cycle
          
        end if

        ! Clean it
        deallocate(R)
        nullify(R)
        
      end do

      deallocate( this%state(is)%orb )
      nullify( this%state(is)%orb )
      
    end do

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

    integer :: io, ind, isite, is, iio
    real(dp) :: occ

    if ( .not. DM%initialized() ) then
      
      ! This routine call will create an initial DM matrix with
      ! minimal sparse pattern. It should only be used
      ! for cases when the atomic DM is the only thing to worry
      ! about. I.e. subsequent usage of the sparse matrix is
      ! ill-adviced for anything but initial Mulliken, initial
      ! density on grids etc.
      call init()
      
    end if

    ! Initialize everything to 0
    DM%M(:) = 0._dp
    
    ! Fill the density matrix with 1's
    do io = 1, this%n_orbital

      isite = this%orbital_site(io)
      is = this%site_state_idx(isite)
      iio = io - this%site_orbital_start(isite) + 1
      occ = this%state(is)%orb(iio)%occ

      ! Figure out if this basis function has initial density
      do ind = sp%rptr(io) , sp%rptr(io) + sp%nrow(io) - 1

        ! Check for the diagonal part
        if ( sp%column(ind) == io ) then
          
          ! Ok, we have a diagonal entry.
          DM%M(ind) = occ
          
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
      call sp%init(this%n_orbital, this%n_orbital, np=1)

      ! Add elements to the sparse matrix
      do io = 1, this%n_orbital
        call sp%add(io, io)
      end do

      ! Finalize the sparse object to remove all unnecessary elements
      call sp%finalize()

      ! Note that now the DM object contains the *only* reference
      ! to sp that exists.
      call DM%init(sp)

    end subroutine init
    
  end subroutine atomic_density_matrix


  !< Calculate the basis function at a position `r` from the center of the basis function
  !<
  !< @param this the basis type
  !< @param specie integer index of the specie
  !< @param r the vector (r_p - r_C) where r_p is the point of evaluation and r_C is the center of the basis functions
  !< @param psi the values of the wave function (has to be large enough to get all)
  subroutine get_psi_all(this, state, r, psi)
    class(basis_ac_t), intent(in) :: this
    integer, intent(in) :: state
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: psi(:)

    ! Local variables
    integer :: norb, io ! number of orbitals, loop
    real(dp) :: dist, FR
    type(orbital_ac_t), pointer :: orb => null()

    norb = this%state(state)%n_orbital
    if ( size(psi) < norb ) then
      stop 'error in get_psi_all'
    end if

    ! Loop orbitals on this state
    dist = sqrt( sum( r ** 2 ) )
    do io = 1, norb

      if ( associated(orb, this%state(state)%orb(io)) ) then
        ! do nothing, the radial value is the same
      else
        orb => this%state(state)%orb(io)
        FR = pspiof_meshfunc_eval(orb%R, dist)
      end if
      
      call grylmr(r(1), r(2), r(3), orb%l, orb%m, psi(io))
      
      psi(io) = psi(io) * FR
      
    end do
    
  end subroutine get_psi_all

  subroutine get_psi_single(this, state, io, r, psi)
    class(basis_ac_t), intent(in) :: this
    integer, intent(in) :: state, io
    real(dp), intent(in) :: r(3)
    real(dp), intent(out) :: psi

    ! Local variables
    real(dp) :: dist
    type(orbital_ac_t), pointer :: orb

    orb => this%state(state)%orb(io)
    call grylmr(r(1), r(2), r(3), orb%l, orb%m, psi)
    dist = sqrt( sum( r ** 2 ) )
    psi = psi * pspiof_meshfunc_eval(orb%R, dist)
    
  end subroutine get_psi_single
  
  !Summary
  !----------------------------------------------------
  subroutine summary(this)
    class(basis_ac_t), intent(in) :: this

    integer :: isite, is, no, io
    character(len=10) :: str
    type(orbital_ac_t), pointer :: orb
    integer, allocatable :: i1(:)
    real(dp), allocatable :: d1(:)

    call yaml_mapping_open("basis_ac")
    
    call yaml_map("Number of sites", this%n_site)

    call yaml_map("Sites first orbital", this%site_orbital_start(:this%n_site))
    call yaml_map("Sites last orbital", this%site_orbital_start(2:this%n_site+1) - 1)

    call yaml_sequence_open("Sites")
    do isite = 1, this%n_site
      
      write(str, '(i0)') isite
      call yaml_sequence(advance="no")
      call yaml_mapping_open(trim(str))

      call yaml_map("Coordinate", this%xyz(:,isite))
      is = this%site_state_idx(isite)
      call yaml_map("State Index", is)
      no = this%state(is)%n_orbital
      call yaml_map("Orbitals", no)

      allocate(i1(no))
      do io = 1, no
        i1(io) = this%state(is)%orb(io)%l
      end do
      call yaml_map("l", i1)
      do io = 1, no
        i1(io) = this%state(is)%orb(io)%m
      end do
      call yaml_map("m", i1)
      deallocate(i1)

      allocate(d1(no))
      do io = 1, no
        d1(io) = this%state(is)%orb(io)%occ
      end do
      call yaml_map("Q0", d1)
      deallocate(d1)
      
      call yaml_mapping_close()
    end do
    call yaml_sequence_close()
    
    call yaml_mapping_close()

  end subroutine summary

end module esl_basis_ac_m
