module mixing_esl
  use prec, only : dp,ip

 implicit none
 private

 public ::               &
            mixing_t,    &
            mixing_linear

 !Data structure for the mixer
 type mixing_t
   real(kind=dp) :: alpha !< Mixing parameter
   contains
    private
    procedure, public :: init
    final  :: cleanup
 end type mixing_t

 contains

   !Initialize the mixer
   !----------------------------------------------------
   subroutine init(this)
     class(mixing_t) :: this
 
     !For the moment this is hardcoded
     this%alpha = 0.3d0

   end subroutine init


   !Release the mixer
   !----------------------------------------------------
   subroutine cleanup(this)
     type(mixing_t) :: this

   end subroutine cleanup

   !Mix the density
   !----------------------------------------------------
   subroutine mixing_linear(this, np, rhoin, rhoout, rhonew)
     type(mixing_t),    intent(in) :: this
     integer,           intent(in) :: np
     real(kind=dp),     intent(in) :: rhoin(:)
     real(kind=dp),     intent(in) :: rhoout(:)
     real(kind=dp),    intent(out) :: rhonew(:)

     integer :: ip

     do ip=1, np
       rhonew(ip) = rhoin(ip)*(1.0d0-this%alpha) + rhoout(ip)*this%alpha
     end do

   end subroutine mixing_linear


end module mixing_esl
