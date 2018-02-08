!< The simplest easy cyclic distribution.
!<
!< The distribution is consecutive and each node gets
!< 1 element in a cyclic manner.
module mpi_dist_cyclic

  use mpi
  use mpi_dist

  implicit none

  private

  ! Possibly we should reduce 
  public :: mpi_dist_cyclic_t

  ! MPI relies on standard integer sizes
  ! However, if one wishes to build this application
  ! with default long integers, we have to be careful.
  integer, parameter :: im_ = selected_int_kind(9)

  ! Data-type for the contained elements
  integer, parameter :: ii_ = selected_int_kind(9)
  
  integer(ii_), parameter :: ONE = 1_ii_
  
  type, extends(mpi_dist_t) :: mpi_dist_cyclic_t
     
   contains
     
     !< Initialize a new MPI distribution
     procedure, public :: init => init_
     
     !< Delete this object (it does NOT disconnect the communicator)
     procedure, public :: delete => delete_
     
     procedure, public :: local_N => local_N_
     procedure, public :: glob_2_loc => global_2_local_
     procedure, public :: loc_2_glob => local_2_global_
     procedure, public :: glob_2_rank => global_2_rank_
     
  end type mpi_dist_cyclic_t
  
contains
  
  subroutine init_(this, comm, global_N)
    class(mpi_dist_cyclic_t), intent(inout) :: this
    !< MPI communicator
    integer(im_), intent(in) :: comm 
    !< Global number of elements
    integer(ii_), intent(in) :: global_N

    call this%mpi_dist_t%init_(comm, global_N)

  end subroutine init_

  subroutine delete_(this)
    class(mpi_dist_cyclic_t), intent(inout) :: this

    call this%mpi_dist_t%delete_()

  end subroutine delete_

  function local_N_(this) result(N)
    class(mpi_dist_cyclic_t), intent(in) :: this
    integer(ii_) :: N

    N = this%global_N / this%size
    ! Figure out if this node has any remaining elements
    if ( this%glob_2_rank(this%global_N) <= this%rank ) then
       N = N + 1
    end if
    
  end function local_N_

  function global_2_local_(this, global) result(local)
    class(mpi_dist_cyclic_t), intent(in) :: this
    integer(ii_), intent(in) :: global
    integer(ii_) :: local
    
    ! Figure out if this node has this element
    if ( this%glob_2_rank(global) == this%rank ) then
       local = global / this%size
    else
       local = -ONE
    end if
  end function global_2_local_

  function local_2_global_(this, local) result(global)
    class(mpi_dist_cyclic_t), intent(in) :: this
    integer(ii_), intent(in) :: local
    integer(ii_) :: global
    global = this%size * local + this%rank
  end function local_2_global_

  function global_2_rank_(this, global) result(rank)
    class(mpi_dist_cyclic_t), intent(in) :: this
    integer(ii_), intent(in) :: global
    integer(ii_) :: rank
    rank = mod(global - ONE, this%global_N)
  end function global_2_rank_

end module mpi_dist_cyclic
