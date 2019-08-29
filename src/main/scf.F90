module esl_scf_m

  use prec, only : dp,ip
  use fdf, only : fdf_get
  
#ifdef WITH_FLOOK
  use esl_flook_global_m, only: LUA
  use dictionary
  use esl_dict_m
  use esl_flook_if_m
#endif

  use esl_basis_m
  use esl_density_m
  use esl_grid_m
  use esl_hamiltonian_m
  use esl_mixing_m
  use esl_potential_m
  use esl_states_m
  use esl_system_m

  ! AC-related functions
  use esl_sparse_pattern_m, only: sparse_pattern_t
  use esl_create_sparse_pattern_ac_m, only: create_sparse_pattern_ac_create
  use esl_overlap_matrix_ac_m, only: overlap_matrix_ac_calculate
  use esl_density_matrix_ac_m, only: density_matrix_ac_next

  implicit none
  private

  public :: scf_t

  !Data structure containing the data for the SCF
  type scf_t
    
    real(dp) :: tol_reldens !< Tolerance for SCF convergence
    integer(ip) :: max_iter !< Maximum number of iterations

    !< Mixer used during the SCF cycle
    type(mixing_t) :: mixer

    !< Density in/out
    type(density_t) :: rho_in, rho_out

    !< Hamiltonian/Potential in/out
    type(hamiltonian_t) :: H

  contains
    private
    procedure, public :: init
    procedure, public :: mix
    
    procedure, public :: loop

    final :: cleanup
  end type scf_t

contains


  !< Initialize density and wfn
  !<
  !< Atomic-Centered orbitals does the following things:
  !< 1. Initialize the sparse pattern by determining couplings based
  !<    purely on the inter-orbital ranges. There will only be
  !<    matrix elements if the distance between two atomic centered
  !<    orbitals are within the sum of their cutoff radius.
  !< 2. Once the sparse pattern has been created we calculate the overlap
  !<    matrix which is fixed for the entirety of the SCF loop.
  !< 3. Finally the Hamiltonian is initialized (an initialization
  !<    of the Hamiltonian also creates the non-SCF matrix elements)
  subroutine init(this, system, states)
    class(scf_t) :: this 
    type(system_t), intent(inout) :: system
    type(states_t), intent(in) :: states

    this%tol_reldens = fdf_get('SCFTolerance',1.0e-6_dp)
    this%max_iter    = fdf_get('SCFMaxIterations', 100)

    call this%mixer%init()

    select case (system%basis%type)
    case ( PLANEWAVES )
      ! TODO anything specific to SCF initialization
    case( ATOMCENTERED )
      ! TODO anything specific to SCF initialization
    end select

    ! This is necessary for both AC and PW
    ! AC only needs the potential class (TODO consider moving the potential_t somewhere else!)
    call this%H%init(system%basis, system%geo, states, &
        periodic=.true.)

    ! Finally we can initialize the densities
    ! This *has* to be done in the end because the density matrices
    ! for AC requires the sparse pattern to be updated.!
    call this%rho_in%init(system)
    call this%rho_out%init(system)

#ifdef WITH_FLOOK
    call esl_dict_var_add('SCF.Mixer.Alpha', this%mixer%alpha)
    call esl_dict_var_add('SCF.Iter.Max', this%max_iter)
    call esl_dict_var_add('SCF.Tolerance', this%tol_reldens)
#endif

  end subroutine init

  !Cleaning up
  !----------------------------------------------------
  subroutine cleanup(this)
    type(scf_t) :: this

  end subroutine cleanup

  !Perform the self-consistent field calculation
  !----------------------------------------------------
  subroutine loop(this, elsi, system, states, smear)
    use yaml_output
    use esl_smear_m
    use esl_states_m
    use esl_elsi_m

    class(scf_t), intent(inout) :: this
    type(elsi_t), intent(inout) :: elsi
    type(system_t), intent(inout) :: system
    type(states_t), intent(inout) :: states
    type(smear_t), intent(inout) :: smear

    integer :: iter !< Interation
    real(dp) :: res

#ifdef WITH_FLOOK
    call esl_dict_var_add('SCF.Residue', res)
#endif

    ! Perform initial guess on the density
    call this%rho_in%guess(system)
    
    !Calc. potentials from guess density
    select case ( system%basis%type )
    case ( PLANEWAVES )
      call this%H%potential%calculate(this%rho_in%pw%density, system%energy)
    case( ATOMCENTERED )
      ! Initialize H0
      call this%H%ac%calculate_H0(system%basis%ac, system%geo)
    end select

    select case (system%basis%type)
    case ( PLANEWAVES )
      ! Randomize the states
      call states%randomize()
    end select

    call yaml_mapping_open("SCF cycle")

    loop_scf: do iter = 1, this%max_iter
      call yaml_map("Iteration", iter)

      ! Since the SCF cycles are *very* different we offload into
      ! separate routines to clarify differences.
      select case ( system%basis%type )
      case ( ATOMCENTERED )
        call SCF_AC()
      case ( PLANEWAVES )
        call SCF_PW()
      end select
      
      ! Calc. energies
      call system%energy%calculate()

      !Test tolerance and print status
      !We use rhonew to compute the relative density
      ! TODO fix the pw/ac part due to the number of states. They are behaving
      !      differently, so perhaps it is better to pass everything?
      !      For now I don't do this...
      res = this%rho_in%residue(system%basis, this%rho_out, states)

      ! Follow convergence
      call yaml_mapping_open("iSCF")
      call yaml_map("Total Energy", system%energy%total)
      call yaml_map("Fermi level", system%energy%fermi)
      call yaml_map("Residue", res)
      call yaml_mapping_close()

#ifdef WITH_FLOOK
      ! Exchange data with Lua, so
      ! Call lua just before we are to initialize a new step.
      call flook_if_call(LUA, LUA_SCF_LOOP)
#endif

      if ( res <= this%tol_reldens ) then
        call yaml_comment("SCF cycle converged.")
        
        exit loop_scf
        
      end if

      ! Perform mixing (in/out)
      call this%mix(system%basis, this%rho_in, this%rho_out)

    end do loop_scf
    
    call yaml_mapping_close()

    call system%energy%display()

  contains

    subroutine SCF_AC()

      real(dp), allocatable :: Vscf(:)

      ! First we calculate the energies for what we know
      system%energy%kinetic = sum(this%H%ac%kin%M(:) * this%rho_in%ac%DM%M(:))
      system%energy%KB = sum(this%H%ac%vkb%M(:) * this%rho_in%ac%DM%M(:))
      
      ! Setup the Hamiltonian so we can solve the eigenvalue problem
      ! First we initialize H0
      ! This takes the pre-calculated *constant* Hamiltonian contributions
      ! and adds them to the SCF Hamiltonian.
      call this%H%ac%setup_H0()
      
      ! Now calculate the charge-density on the grid to be able to
      ! calculate contributions from the potentials to the SCF Hamiltonian
      ! This assumes DM stored in rho_in is the latest available!
      call this%rho_in%calculate(system)

      ! Calculate potentials and then add to H
      ! Now call the potentials_t%calculate which does:
      !   1. Calculate the XC potential.
      !   2. Calculate the Hartree potential
      !   3. Calculate the external local potential
      call this%H%potential%calculate(this%rho_in%ac%rho, system%energy)

      ! Prepare the Vscf to add potential to H
      allocate(Vscf(system%basis%ac%grid%np))
      Vscf(:) = this%H%potential%Hartree(:) + this%H%potential%Vxc(:) + this%H%potential%external(:)
      call this%H%ac%add_potential(system%basis%ac, Vscf)
      deallocate(Vscf)

      ! Diagonalization (ELSI/KSsolver)
      call this%H%eigensolver(system%basis, states)

      ! Update occupations
      call smear%calc_fermi_occ(elsi, states)
      ! Copy over fermi-level
      system%energy%fermi = smear%fermi_level
      
      ! Calculate new density matrix from states
      ! rho%calculate has system and states as optional
      ! This allows the same routine to calculate
      ! 1. the DM expanded to the grid (if system present)
      ! 2. the DM from the states (if system not present, and states present)
      !
      ! TODO consider moving this to a more clear implementation
      ! I.e. currently rho%calculate can do both grid and DM calculations
      ! for AC. However, doing it like this ensures that we only have one method.... :(
      call this%rho_out%calculate(system, states)

      ! H is now the input and DM is DM(H)
      ! Thus we can calculate the eigenvalue energy now
      !TODO for spin
      system%energy%eigenvalues = sum(this%H%ac%H(1)%M(:) * this%rho_out%ac%DM%M(:))

    end subroutine SCF_AC

    subroutine SCF_PW()
      ! Diagonalization (ELSI/KSsolver)
      call this%H%eigensolver(system%basis, states)
      
      ! Update occupations
      call smear%calc_fermi_occ(elsi, states)
      ! Copy over fermi-level
      system%energy%fermi = smear%fermi_level
      
      ! Calculate density from the new states
      call this%rho_out%calculate(system, states)
      
      ! Calculate necessary potentials
      call this%H%potential%calculate(this%rho_out%pw%density, system%energy)

      ! TODO Check that the energies are correct. Added eigenvalues times occupations!
      system%energy%eigenvalues = sum(states%eigenvalues * states%occ_numbers)

    end subroutine SCF_PW

  end subroutine loop
  
  ! Perform mixing step
  !----------------------------------------------------
  subroutine mix(this, basis, in, out)
    class(scf_t) :: this
    type(basis_t),  intent(in)   :: basis
    type(density_t), intent(inout) :: in, out 

    integer :: np

    select case ( basis%type )
    case ( PLANEWAVES )

      np = basis%pw%grid%np
      call this%mixer%linear(np, in%pw%density(1:np), &
          out%pw%density(1:np), in%pw%density(1:np))
      
    case ( ATOMCENTERED )

      ! Since the linear mixing does not do look-ahead calls we can alias
      ! the same array
      np = in%ac%DM%sp%nz
      call this%mixer%linear(np, in%ac%DM%M, out%ac%DM%M, in%ac%DM%M)
      
    end select

  end subroutine mix

end module esl_scf_m
