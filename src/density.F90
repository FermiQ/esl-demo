module density_esl
  use prec, only : dp,ip

  use basis_esl
  use grid_esl

  implicit none
  private

  public :: density_t
  
  !Data structure for the density
  type density_t
    integer :: np !< Copied from grid

    real(kind=dp), allocatable :: density(:)
    real(kind=dp), allocatable :: density_matrix(:,:) 
   contains
     procedure, public :: init
     procedure, public :: calculate
     procedure, public :: get_den
     procedure, public :: set_den
     final  :: cleanup
  end type density_t

  contains

   !Initialize the density
   !----------------------------------------------------
   subroutine init(this, basis, grid)
     class(density_t), intent(inout) :: this
     type(basis_t),   intent(in)    :: basis
     type(grid_t),    intent(in)    :: grid

     !Parse the informations from the input file
     select case(basis%basis_type)
       case(PLANEWAVES)
         !Initialization structures for the PW case
         allocate(this%density(1:grid%np))
         this%density(1:grid%np) = 0.d0
       case(ATOMICORBS)
         !Initialization structures for the LO case
         !TEMP
         allocate(this%density(1:grid%np))
         this%density(1:grid%np) = 0.d0
     end select 

     this%np = grid%np

   end subroutine init
 
   !Release
   !----------------------------------------------------
   subroutine cleanup(this)
     type(density_t), intent(inout) :: this

     if(allocated(this%density_matrix)) deallocate(this%density_matrix)
     if(allocated(this%density)) deallocate(this%density)

   end subroutine cleanup

   !Calc density
   !----------------------------------------------------
   subroutine calculate(this, basis)
     class(density_t), intent(inout) :: this
     type(basis_t),   intent(in)    :: basis

     select case(basis%basis_type)
       case(PLANEWAVES)

       case(ATOMICORBS)

     end select

   end subroutine calculate


   !Copy the density to an array
   !----------------------------------------------------
   subroutine get_den(this, rho)
     class(density_t) :: this
     real(kind=dp),    intent(out) :: rho(:)

     integer :: ip

     forall(ip=1:this%np)
       rho(ip) = this%density(ip)
     end forall

   end subroutine get_den

   !Copyi the density from an array
   !----------------------------------------------------
   subroutine set_den(this, rho)
     class(density_t) :: this
     real(kind=dp),    intent(in) :: rho(:)

     integer :: ip

     forall(ip=1:this%np)
       this%density(ip) = rho(ip)
     end forall

   end subroutine set_den


end module density_esl
