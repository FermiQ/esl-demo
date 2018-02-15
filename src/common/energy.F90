module esl_energy_m
  use prec, only : dp

  use esl_states_m

  implicit none
  private

  public ::                          &
      energy_t

  !Data structure for the energy
  type energy_t
    real(dp) :: total
    real(dp) :: eigenvalues
    real(dp) :: hartree
    real(dp) :: ionion
    real(dp) :: extern
    real(dp) :: exchange
    real(dp) :: correlation
    real(dp) :: int_nvxc
    real(dp) :: kinetic
    real(dp) :: entropy
    real(dp) :: fermi

  contains
    private
    procedure, public :: init
    procedure, public :: calculate
    procedure, public :: display
  end type energy_t


contains

  !Initialize the energies
  !----------------------------------------------------
  subroutine init(this)
    class(energy_t) :: this

    this%total = 0._dp
    this%eigenvalues = 0._dp
    this%hartree = 0._dp
    this%ionion = 0._dp
    this%extern = 0._dp
    this%exchange = 0._dp   
    this%correlation = 0._dp
    this%int_nvxc = 0._dp
    this%kinetic = 0._dp
    this%entropy = 0._dp
    this%fermi = 0._dp

  end subroutine init


  !Compute the total energy
  !----------------------------------------------------
  subroutine calculate(this, states)
    class(energy_t) :: this
    type(states_t), intent(in) :: states

    this%eigenvalues=sum(states%eigenvalues)

    this%total = this%ionion +  this%eigenvalues  + this%extern &
                 -this%hartree + this%exchange + this%correlation + this%kinetic - this%int_nvxc

  end subroutine calculate


  !Display de different components of the total energy
  !----------------------------------------------------
  subroutine display(this)
    use yaml_output
    class(energy_t) :: this

    call yaml_mapping_open("Energy", advance='NO')
    call yaml_comment("Hartree", hfill = "-")
    call yaml_map("Total", this%total)
    call yaml_map("Eigenvalues", this%eigenvalues)
    call yaml_map("Hartree", this%hartree)
    call yaml_map("Ion-ion", this%ionion)
    call yaml_map("Extern", this%extern)
    call yaml_map("Int_nvxc", this%int_nvxc)
    call yaml_map("Kinetic", this%kinetic)
    call yaml_map("Fermi-level", this%fermi)
    call yaml_mapping_close()

  end subroutine display

end module esl_energy_m
