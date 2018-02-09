module forces_esl
  use prec, only : dp,ip

  implicit none
  private

  public ::                          &
            forces_t
  
  !Data structure for the energy
  type forces_t
    real(kind=dp), allocatable :: total(:,:)
    real(kind=dp), allocatable :: loc(:,:)
    real(kind=dp), allocatable :: nl(:,:)
    real(kind=dp), allocatable :: ionion(:,:)
    contains
    private
    procedure, public :: init
    final :: cleanup
  end type forces_t


  contains

   !Initialize the forcess
   !----------------------------------------------------
   subroutine init(this, natoms)
     class(forces_t) :: this

     allocate(this%total(1:3,1:natoms))
     this%total(1:3,1:natoms) = 0.d0
     allocate(this%loc(1:3,1:natoms))
     this%loc(1:3,1:natoms) = 0.d0
     allocate(this%nl(1:3,1:natoms))
     this%nl(1:3,1:natoms) = 0.d0
     allocate(this%ionion(1:3,1:natoms))
     this%ionion(1:3,1:natoms) = 0.d0

   end subroutine init

   !Release
   !----------------------------------------------------
   subroutine cleanup(this)
     type(density_t), intent(inout) :: this

     if(allocated(this%total)) deallocate(this%total)
     if(allocated(this%loc)) deallocate(this%loc)
     if(allocated(this%nl)) deallocate(this%nl)
     if(allocated(this%ionion)) deallocate(this%ionion)
   end subroutine cleanup

 
end module forces_esl
