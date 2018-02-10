module esl_energy_t
  use prec, only : dp

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

  end subroutine init


  !Compute the total energy
  !----------------------------------------------------
  subroutine calculate(this)
    class(energy_t) :: this

    this%total = this%ionion + this%correlation 
    this%total = this%total + 0.5_dp*(this%eigenvalues + this%kinetic + this%extern + this%int_nvxc)

  end subroutine calculate


  !Display de different components of the total energy
  !----------------------------------------------------
  subroutine display(this)
    use yaml_output
    class(energy_t) :: this

    integer :: i

    call yaml_mapping_open("Energy")
    call yaml_map("Total", this%total)
    call yaml_map("Eigenvalues", this%eigenvalues)
    call yaml_map("Hatree", this%hartree)
    call yaml_map("Ion-ion", this%ionion)
    call yaml_map("Extern", this%extern)
    call yaml_map("Int_nvxc", this%int_nvxc)
    call yaml_map("Kinetic", this%kinetic)
    call yaml_mapping_close()

  end subroutine display


end module esl_energy_t
