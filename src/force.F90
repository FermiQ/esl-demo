module esl_force_t

  use prec, only : dp

  implicit none
  private

  public :: force_t

  !Data structure for the energy
  type force_t
     
     real(dp), allocatable :: total(:,:)
     real(dp), allocatable :: loc(:,:)
     real(dp), allocatable :: nl(:,:)
     real(dp), allocatable :: ionion(:,:)
     
   contains
     private
     
     procedure, public :: init
     final :: cleanup
     
  end type force_t

contains

  !Initialize the forces
  !----------------------------------------------------
  subroutine init(this, natoms)
    class(force_t), intent(inout) :: this
    integer, intent(in) :: natoms

    allocate(this%total(1:3,1:natoms))
    this%total(1:3,1:natoms) = 0._dp
    allocate(this%loc(1:3,1:natoms))
    this%loc(1:3,1:natoms) = 0._dp
    allocate(this%nl(1:3,1:natoms))
    this%nl(1:3,1:natoms) = 0._dp
    allocate(this%ionion(1:3,1:natoms))
    this%ionion(1:3,1:natoms) = 0._dp

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(force_t), intent(inout) :: this

    if(allocated(this%total)) deallocate(this%total)
    if(allocated(this%loc)) deallocate(this%loc)
    if(allocated(this%nl)) deallocate(this%nl)
    if(allocated(this%ionion)) deallocate(this%ionion)
    
  end subroutine cleanup

end module esl_force_t
