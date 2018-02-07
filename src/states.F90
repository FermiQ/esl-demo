module states_esl
  use prec, only : dp,ip

  use basis_esl

 implicit none
 private

 public ::                   &
           states_t,         &
           states_init,      &
           states_end
          
 !Data structure for the states
 type states_t
   integer :: nstates
   integer :: nkpt
   integer :: nspin

   type(wfn_t), allocatable :: states(:,:,:)  !nstates, nspin, nkpt
 end type states_t

 type wfn_t
   real(kind=dp), allocatable :: coef(:) !<Coefficients of the wavefunction in the basis
 end type wfn_t

 contains

   !Initialize the states
   !----------------------------------------------------
   subroutine states_init(this, basis, nstates, npsin, nkpt)
     type(states_t), intent(inout) :: this
     type(basis_t),  intent(in)    :: basis
     integer,        intent(in)    :: nstates
     integer,        intent(in)    :: nspin
     integer,        intent(in)    :: nkpt

     integer :: ist, isp, ik

     this%nstates = nstates
     this%nspin = nspin
     this%nkpt = nkpt

     allocate(this%states(1:nstates, 1:nspin, 1:nkpt))
     
     do ik = 1, nkpt
       do isp = 1, nspin
         do ist = 1, nstates
           allocate(this%states(ist, isp, ik)%coef(1:basis%size)
         end do 
       end do
     end do

   end subroutine states_init


   !Release the states
   !----------------------------------------------------
   subroutine states_end(this)
     type(states_t):: this

     allocate(this%states(1:nstates, 1:nspin, 1:nkpt))

     do ik = 1, nkpt
       do isp = 1, nspin
         do ist = 1, nstates
           if(this%states(ist, isp, ik)%coef) &
             deallocate(this%states(ist, isp, ik)%coef)
         end do 
       end do
     end do

     if(allocated(this%states)) deallocate(this%states)

   end subroutine states_end


end module potential_esl
