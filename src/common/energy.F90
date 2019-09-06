
!< Energy module implementing a type to retain energies.
!<
!< The type should be used for energies calculated.
!< Generally the energies only gets updated upon explicit calls.
!< So interaction with this type requires knowledge of the used
!< energies calculation.
!< The total energy is an important energy that only gets updated
!< upon `call energy%calculate(states)` calls.
!< Thus using the total energy requires that `calculate` has been
!< called.
!<
!< All energies are in Hartree.
module esl_energy_m
  use prec, only : dp

  implicit none
  private

  public :: energy_t

  !< Data structure for energy terms
  type energy_t
    !< Total energy, only updated on `calculate` call
    real(dp) :: total = 0._dp
    !< Energy from KS eigenvalues
    real(dp) :: eigenvalues = 0._dp
    !< Hartree potential energy
    real(dp) :: hartree = 0._dp
    !< Kleynman-Bylander energies
    real(dp) :: KB = 0._dp
    !< Ion-ion interaction energy
    real(dp) :: ionion = 0._dp
    !< External energy
    real(dp) :: extern = 0._dp
    !< Exchange energy
    real(dp) :: exchange = 0._dp
    !< Correlation energy
    real(dp) :: correlation = 0._dp
    !< TODO ?
    real(dp) :: int_nvxc = 0._dp
    !< Kinetic energy
    real(dp) :: kinetic = 0._dp
    !< Entropy (not energy, unit-less)
    real(dp) :: entropy = 0._dp
    !< Fermi-level of system
    real(dp) :: fermi = 0._dp

  contains
    private
    
    procedure, public :: init
    procedure, public :: calculate
    procedure, public :: display
    
  end type energy_t


contains

  !< Initialize data type by setting all contained variables to 0.
  subroutine init(this)
    class(energy_t) :: this

    this%total = 0._dp
    this%eigenvalues = 0._dp
    this%hartree = 0._dp
    this%KB = 0._dp
    this%ionion = 0._dp
    this%extern = 0._dp
    this%exchange = 0._dp   
    this%correlation = 0._dp
    this%int_nvxc = 0._dp
    this%kinetic = 0._dp
    this%entropy = 0._dp
    this%fermi = 0._dp

  end subroutine init


  !< Update total energy
  subroutine calculate(this)
    class(energy_t) :: this

    this%total = this%ionion + this%eigenvalues + this%extern &
                 - this%hartree + this%exchange + this%correlation + this%kinetic - this%int_nvxc

  end subroutine calculate


  !< Print to std-out the energy decomposition
  subroutine display(this)
    use yaml_output
    class(energy_t) :: this
    
    call yaml_mapping_open("Energy", advance='NO')
    call yaml_comment("Hartree", hfill = "-")
    call yaml_map("Total", this%total)
    call yaml_map("Eigenvalues", this%eigenvalues)
    call yaml_map("Hartree", this%hartree)
    call yaml_map("KB", this%KB)
    call yaml_map("Ion-ion", this%ionion)
    call yaml_map("Extern", this%extern)
    call yaml_map("Exchange", this%exchange)
    call yaml_map("Correlation", this%correlation)
    call yaml_map("Int_nvxc", this%int_nvxc)
    call yaml_map("Kinetic", this%kinetic)
    call yaml_map("Fermi-level", this%fermi)
    call yaml_mapping_close()

  end subroutine display

end module esl_energy_m
