module states_esl
  use prec, only : dp,ip

  use basis_esl
  use numeric_esl

 implicit none
 private

 public ::                   &
           states_t,         &
           states_init,      &
           states_end,       &
           states_randomize

 type wfn_t
   real(kind=dp),    allocatable :: dcoef(:) !<Coefficients of the wavefunction in the basis
   complex(kind=dp), allocatable :: zcoef(:)
 end type wfn_t

          
 !Data structure for the states
 type states_t
   integer :: nstates
   integer :: nkpt
   integer :: nspin
   integer :: ncoef

   logical :: complex_states

   type(wfn_t), allocatable :: states(:,:,:)  !nstates, nspin, nkpt
   real(kind=dp), allocatable :: occ_numbers(:,:,:)
 end type states_t

 contains

   !Initialize the states
   !----------------------------------------------------
   subroutine states_init(this, basis, nstates, nspin, nkpt)
     type(states_t), intent(inout) :: this
     type(basis_t),  intent(in)    :: basis
     integer,        intent(in)    :: nstates
     integer,        intent(in)    :: nspin
     integer,        intent(in)    :: nkpt

     integer :: ist, isp, ik

     this%nstates = nstates
     this%nspin = nspin
     this%nkpt = nkpt
     this%ncoef = basis%size

     allocate(this%states(1:nstates, 1:nspin, 1:nkpt))
     allocate(this%occ_numbers(1:nstates, 1:nspin, 1:nkpt))
     this%occ_numbers(1:nstates, 1:nspin, 1:nkpt) = 0._dp

     select case(basis%basis_type)
       case(PLANEWAVES)
         this%complex_states = .true.
         do ik = 1, nkpt
           do isp = 1, nspin
             do ist = 1, nstates
               allocate(this%states(ist, isp, ik)%zcoef(1:basis%size))
             end do 
           end do
         end do
       case(ATOMICORBS)
         this%complex_states = .false.
         do ik = 1, nkpt
           do isp = 1, nspin
             do ist = 1, nstates
               allocate(this%states(ist, isp, ik)%dcoef(1:basis%size))
             end do
           end do
         end do
     end select

   end subroutine states_init


   !Release the states
   !----------------------------------------------------
   subroutine states_end(this)
     type(states_t):: this

     integer :: ist, isp, ik

     if(allocated(this%states)) then
        do ik = 1, this%nkpt
           do isp = 1, this%nspin
              do ist = 1, this%nstates
                 if(allocated(this%states(ist, isp, ik)%dcoef)) &
                      deallocate(this%states(ist, isp, ik)%dcoef)
                 if(allocated(this%states(ist, isp, ik)%zcoef)) &
                      deallocate(this%states(ist, isp, ik)%zcoef)
              end do
           end do
        end do
        deallocate(this%states)
     end if
     if(allocated(this%occ_numbers)) deallocate(this%occ_numbers)

   end subroutine states_end


   !Randomize the states
   !----------------------------------------------------
   subroutine states_randomize(this)
     type(states_t):: this

     integer :: ist, isp, ik
     real(kind=dp), allocatable :: tmp_re(:), tmp_im(:)

     call init_random()

     if(this%complex_states) then
       allocate(tmp_re(1:this%ncoef))
       allocate(tmp_im(1:this%ncoef))
       do ik = 1, this%nkpt
         do isp = 1, this%nspin
           do ist = 1, this%nstates
             call random_number(tmp_re(1:this%ncoef))
             call random_number(tmp_im(1:this%ncoef))
             this%states(ist, isp, ik)%zcoef(1:this%ncoef) = tmp_re(1:this%ncoef) + cmplx(0.d0,1.d0)*tmp_im(1:this%ncoef) 
           end do
         end do
       end do
       deallocate(tmp_re, tmp_im)
     else 
       do ik = 1, this%nkpt
         do isp = 1, this%nspin
           do ist = 1, this%nstates
             call random_number(this%states(ist, isp, ik)%dcoef(1:this%ncoef))
           end do
         end do
       end do
     end if

   end subroutine states_randomize

end module states_esl
