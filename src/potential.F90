module potential_esl
  use prec, only : dp,ip

  use basis_esl
  use density_esl
  use energy_esl
  use psolver_esl
  use states_esl

 implicit none
 private

 public ::                   &
           potential_t,      &
           potential_init,   &
           potential_end,    &
           potential_calc
          
 !Data structure for the potentials
 type potential_t
   integer :: np !< Number of points in real space

   real(kind=dp), allocatable :: hartree(:)  !Hartree potential
   real(kind=dp), allocatable :: external(:) !External local potential
   real(kind=dp), allocatable :: xc(:,:)  !xc potential

   real(dp) :: ionicOffset !< Offset of the external potential

   type(psolver_t) :: psolver
 end type 

 contains

   !Initialize the potentials
   !----------------------------------------------------
   subroutine potential_init(this, basis, states)
     type(potential_t) :: this
     type(basis_t), intent(in) :: basis
     type(states_t), intent(in):: states

     integer :: ndim
     real(kind=dp) :: hgrid
     character(len = 1) :: geocode

     this%np = ndim*ndim*ndim
     

     allocate(this%hartree(1:this%np))
     allocate(this%external(1:this%np))
     allocate(this%xc(1:this%np, 1:states%nspin))

     select case(basis%basis_type)
       case(PLANEWAVES)
         geocode = 'P'
       case(ATOMICORBS)
         geocode = 'F'
      end select

      this%ionicOffset = 0._dp
     call this%psolver%init(1, 1, geocode, (/ndim, ndim, ndim/), (/hgrid, hgrid, hgrid/))    

     !Here we need to init the libxc and pspio parts 

   end subroutine potential_init


   !Release the potentials
   !----------------------------------------------------
   subroutine potential_end(this)
     type(potential_t):: this

     if(allocated(this%hartree)) call deallocate(this%hartree)
     if(allocated(this%external)) call deallocate(this%external)
     if(allocated(this%xc)) call deallocate(this%xc)

   end subroutine potential_end


   !Compute the different potentials
   !----------------------------------------------------
   subroutine potential_calc(this, density, energy)
     type(potential_t), intent(inout) :: this
     type(density_t),   intent(in)    :: density
     type(energy_t),    intent(inout) :: energy

     call this%psolver%h_potential(density, this%hartree, this%np, &
          & this%external, this%ionicOffset, energy%hartree)

     !Here we need to compute the xc potential

   end subroutine potential_calc

end module potential_esl
