module esl_scf_m

  use prec, only : dp,ip
  use fdf, only : fdf_get

  use esl_basis_m
  use esl_density_m
  use esl_grid_m
  use esl_hamiltonian_m
  use esl_mixing_m
  use esl_potential_m
  use esl_system_m

  implicit none
  private

  public :: &
      scf_t

  !Data structure containing the data for the SCF
  type scf_t 
    real(dp) :: tol_reldens 
    integer(ip) :: max_iter !< Maximum number of iterations

    !< Mixer used during the SCF cycle
    type(mixing_t) :: mixer

    !< Density in/out
    type(density_t) :: rho_in, rho_out

    !< Hamiltonian/Potential in/out
    type(hamiltonian_t) :: H_in, H_out

  contains
    private
    procedure, public :: init
    procedure, public :: mix
    
    procedure, public :: loop

    final :: cleanup
  end type scf_t

contains


  !Initialize density and wfn
  !----------------------------------------------------
  subroutine init(this)
    class(scf_t)  :: this

    this%tol_reldens = fdf_get('SCFTolerance',1.0e-6_dp)
    this%max_iter    = fdf_get('SCFMaxIterations', 100)
    !Parse here the data for the SCF

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
    use elsi_wrapper_esl

    class(scf_t),  intent(inout) :: this
    type(elsi_t),  intent(inout) :: elsi
    type(system_t),   intent(in) :: system
    type(states_t),intent(inout) :: states
    type(smear_t), intent(inout) :: smear

    integer :: iter !< Interation
    real(dp) :: res

    ! Initialize the densities
    call this%rho_in%init(system%basis)
    call this%rho_out%init(system%basis)

    ! Perform initial guess on the density
    call this%rho_in%guess(system%basis, system%geo)

    call yaml_mapping_open("SCF cycle")

    loop_scf: do iter = 1, this%max_iter
      call yaml_map("Iteration", iter)

      ! Diagonalization (ELSI/KSsolver)
      call this%H_in%eigensolver(system%basis, states)

      ! Update occupations
      call smear_calc_fermi_and_occ(smear, elsi, states)

      ! Calculate density
      call this%rho_in%calculate(system%basis, states, out=this%rho_out)
      
      !Calc. potentials
!      call this%H_in%potentials%calculate(this%rho_out, this%H_in%energy)

      !Calc. energies
      call this%H_in%energy%calculate()  

      !Test tolerance and print status
      !We use rhonew to compute the relative density
      ! TODO fix the pw/ac part due to the number of states. They are behaving
      !      differently, so perhaps it is better to pass everything?
      !      For now I don't do this...
      res = this%rho_in%residue(system%basis, this%rho_out, states)
      call yaml_map("Residue", res)
!      call yaml_map("Rel. Density", reldens)
      if ( res <= this%tol_reldens ) then
        
        call yaml_comment("SCF cycle converged.")
        
        exit loop_scf
        
      end if

      ! Perform mixing
      call this%mix(system%basis, this%rho_in, this%rho_out)

      !Update Hamiltonian matrix

    end do loop_scf
    
    call yaml_mapping_close()

    call this%H_out%energy%display()

  end subroutine loop
  
  ! Perform mixing step
  !----------------------------------------------------
  subroutine mix(this, basis, rho_in, rho_out)
    class(scf_t) :: this
    type(basis_t),  intent(in)   :: basis
    type(density_t), intent(out) :: rho_in, rho_out 

    real(kind=dp), allocatable :: rhonew(:)

    select case (basis%type)
    case ( PLANEWAVES )  
      allocate(rhonew(1:rho_in%np))
      call this%mixer%linear(rho_in%np, rho_in%density_pw%density, rho_out%density_pw%density, rhonew)
      call rho_in%density_pw%set_den(rhonew)
      deallocate(rhonew)
    case ( ATOMCENTERED )
    !TODO
!!$      call mixing_linear(this%mixer, this%np, this%rhoin, this%rhoout, this%rhonew)
!!$      call this%ac%set_den(this%rhonew)
    end select

  end subroutine mix

end module esl_scf_m
