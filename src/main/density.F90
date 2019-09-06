module esl_density_m

  use prec, only : dp
  use esl_basis_m
  use esl_elsi_m
  use esl_density_ac_m
  use esl_density_pw_m
  use esl_geometry_m
  use esl_mixing_m
  use esl_hamiltonian_m
  use esl_states_m
  use esl_system_m

  implicit none

  private

  public :: density_t

  !Data structure for the density
  type density_t
    type(density_ac_t) :: ac
    type(density_pw_t) :: pw
    
  contains
    
    private
    procedure, public :: init
    procedure, public :: guess
    procedure, public :: calculate
    procedure, public :: calculate_density_matrix
    procedure, public :: residue
    final  :: cleanup
    
  end type density_t

contains

  !< Initialize density objects
  subroutine init(this, system)
    class(density_t), intent(inout) :: this
    type(system_t), intent(in) :: system

    select case ( system%basis%type )
    case ( PLANEWAVES )
      
      call this%pw%init(system%basis%pw)
      
    case ( ATOMCENTERED )
      
      call this%ac%init(system%basis%ac)
      
    end select

  end subroutine init

  !Guess the initial density from the atomic orbitals
  !----------------------------------------------------
  subroutine guess(this, system)
    class(density_t), intent(inout) :: this
    type(system_t), intent(in) :: system
   
    select case ( system%basis%type )
    case ( PLANEWAVES )

      ! TODO fix interface 
      call this%pw%guess(system%geo)
      
    case ( ATOMCENTERED )

      ! Guess the initial density from the atomic fillings
      call this%ac%guess(system%basis%ac)
      
    end select

  end subroutine guess

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(density_t), intent(inout) :: this

  end subroutine cleanup

  ! Calculate new density based on a system and states
  !----------------------------------------------------
  subroutine calculate(this, system, states)
    use yaml_output, only: yaml_map
    class(density_t), intent(inout) :: this
    type(system_t), intent(inout) :: system
    type(states_t), intent(in) :: states

    ! Calculate density on the grid
    select case ( system%basis%type )
    case ( PLANEWAVES )

      ! Calculate density
      call this%pw%calculate(states)
      call yaml_map("Norm", system%basis%pw%grid%integrate(this%pw%density))

    case ( ATOMCENTERED )

      ! We are requesting the current DM expanded onto the grid
      call this%ac%calculate(system%basis%ac)
      
    end select

  end subroutine calculate

  ! Calculate new density matrix based on the states
  !----------------------------------------------------
  subroutine calculate_density_matrix(this, system, states)
    use yaml_output, only: yaml_map
    class(density_t), intent(inout) :: this
    type(system_t), intent(inout) :: system
    type(states_t), intent(in) :: states

    ! Calculate density on the grid
    select case ( system%basis%type )
    case ( ATOMCENTERED )

      call this%ac%calculate_density_matrix(states)
      
    end select

  end subroutine calculate_density_matrix


  !< Calculate the relative difference between two densities
  function residue(this, basis, other, states) result(res)
    class(density_t),    intent(in) :: this, other
    type(basis_t),       intent(in) :: basis
    type(states_t),      intent(in) :: states

    real(dp) :: res

        ! Calculate density on the grid
    select case ( basis%type )
    case ( PLANEWAVES )

      res = this%pw%residue(basis%pw%grid, states%nel, other%pw)
      
    case ( ATOMCENTERED )

      res = this%ac%residue(other%ac)

    end select

  end function residue

end module esl_density_m
