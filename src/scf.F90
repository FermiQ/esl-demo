module esl_scf_m

  use prec, only : dp,ip
  use fdf, only : fdf_get

  use basis_esl
  use esl_density_t
  use esl_grid_m
  use esl_hamiltonian_m
  use esl_mixing_m
  use esl_potential_m
  use esl_system_m

  implicit none
  private

  public ::                       &
       scf_t,                &
       scf_loop

  !Data structure containing the data for the SCF
  type scf_t 
     real(kind=dp) :: tol_reldens 
     integer(kind=ip) :: max_iter !< Maximum number of iterations

     type(mixing_t)   :: mixer
   contains
     private
     procedure, public :: init
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

    call this%mixer%init()

  end subroutine init

  !Cleaning up
  !----------------------------------------------------
  subroutine cleanup(this)
    type(scf_t) :: this


  end subroutine cleanup

  !Perform the self-consistent field calculation
  !----------------------------------------------------
  subroutine scf_loop(this, elsi, hamiltonian, system, states, smear)
    use yaml_output
    use esl_smear_m
    use esl_states_m
    use elsi_wrapper_esl

    type(scf_t),         intent(inout) :: this
    type(elsi_t), intent(inout) :: elsi
    type(hamiltonian_t), intent(inout) :: hamiltonian
    type(system_t),         intent(in) :: system
    type(states_t),         intent(in) :: states
    type(smear_t), intent(inout) :: smear

    integer :: iter, ip !< Interation
    real(kind=dp), allocatable :: rhoin(:)
    real(kind=dp), allocatable :: rhoout(:)
    real(kind=dp), allocatable :: rhonew(:)
    real(kind=dp) :: reldens

    allocate(rhoin(1:system%grid%np))
    allocate(rhoout(1:system%grid%np))
    allocate(rhonew(1:system%grid%np))

    call hamiltonian%density%guess(system)

    call yaml_mapping_open("SCF cycle")

    do iter = 1, this%max_iter
       call yaml_map("Iteration", iter)

       !Diagonalization (ELSI/KSsolver)

       !Update occupations
       call smear_calc_fermi_and_occ(smear, elsi, states)

       !Saving the in density for the mixing
       call hamiltonian%density%get_den(rhoin)
       !Calc. density
       call hamiltonian%density%calculate(system%basis)
       !Saving the out density for the mixing
       call hamiltonian%density%get_den(rhoout)
       !Calc. potentials
       call hamiltonian%potentials%calculate(hamiltonian%density, hamiltonian%energy)

       !Calc. energies
       call hamiltonian%energy%calculate()  

       !Test tolerance and print status
       !We use rhonew to compute the relative density
       do ip = 1, system%grid%np
          rhonew(ip) = abs(rhoout(ip) - rhoin(ip))
       end do
       call integrate(system%grid, rhonew, reldens)
       reldens = reldens/real(states%nel)
       call yaml_map("Rel. Density", reldens)
       if(reldens <= this%tol_reldens) then
          call yaml_comment("SCF cycle converged.")
          exit
       end if

       !Mixing (BLAS/LAPACK)
       call mixing_linear(this%mixer, system%grid%np, rhoin, rhoout,rhonew)  
       call hamiltonian%density%set_den(rhonew)

       !Update Hamiltonian matrix

    end do
    call yaml_mapping_close()

    call hamiltonian%energy%display()

    deallocate(rhoin, rhoout, rhonew)

  end subroutine scf_loop

end module esl_scf_m
