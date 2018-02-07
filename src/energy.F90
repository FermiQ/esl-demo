module energy_esl
  use prec, only : dp,ip

  implicit none
  private

  public ::                          &
            energy_t
  
  !Data structure for the energy
  type energy_t
    real(kind=dp) :: total
    real(kind=dp) :: eigenvalues
    real(kind=dp) :: hartree
    real(kind=dp) :: ionion
    real(kind=dp) :: extern
    real(kind=dp) :: exchange
    real(kind=dp) :: correlation
    real(kind=dp) :: int_nvxc
    real(kind=dp) :: kinetic
    real(kind=dp) :: entropy
    contains
    private
    procedure, public :: init
    procedure, public :: calculate
  end type energy_t


  contains

   !Initialize the energies
   !----------------------------------------------------
   subroutine calculate(this)
     class(energy_t) :: this

     this%total = 0.d0
     this%eigenvalues = 0.d0
     this%hartree = 0.d0
     this%ionion = 0.d0
     this%extern = 0.d0
     this%exchange = 0.d0   
     this%correlation = 0.d0
     this%int_nvxc = 0.d0
     this%kinetic = 0.d0
     this%entropy = 0.d0

   end subroutine calculate


   !Compute the total energy
   !----------------------------------------------------
   subroutine calculate(this)
     class(energy_t) :: this

     this%total = this%ionion + this%correlation 
     this%total = this%total + 0.5d0*(this%eigenvalues + this%kinetic + this%extern + this%int_nvxc)

   end subroutine calculate
 
end module energy_esl
