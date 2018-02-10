module xc_esl
  use prec, only : dp

  use esl_density_t


  implicit none
  private

  public :: xc_t
  
  !Data structure for the xc potential
  type xc_t
    integer :: exchange
    integer :: correlation
   contains
     procedure, public :: init
     procedure, public :: calculate
     final  :: cleanup
  end type xc_t

  contains

   !Initialize the xc 
   !----------------------------------------------------
   subroutine init(this)
     class(xc_t), intent(inout) :: this

     !Here parse the exchange and correlation functionals
     !and do libxc related initialization

   end subroutine init
 
   !Release
   !----------------------------------------------------
   subroutine cleanup(this)
     type(xc_t), intent(inout) :: this

   end subroutine cleanup

   !Calc the xc potential from the density
   !----------------------------------------------------
   subroutine calculate(this, density, vxc)
     class(xc_t),           intent(in) :: this
     type(density_t),       intent(in) :: density 
     real(kind=dp),        intent(out) :: vxc(:,:)

     !Here add the libxc business

   end subroutine calculate


end module xc_esl
