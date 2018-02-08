!< The simplest easy cyclic distribution.
!<
!< The distribution is consecutive and each node gets
!< 1 element in a cyclic manner.
module mpi_dist_cyclic

#ifdef WITH_MPI
  use mpi
#endif
  use mpi_dist

  ! MPI relies on standard integer sizes
  ! However, if one wishes to build this application
  ! with default long integers, we have to be careful.
  use prec, only: im_ => ip
  
  ! To increase precision for the container
  ! simply change the pointed to variable.
  use prec, only: ii_ => ip
  
  implicit none

  private

  ! Possibly we should reduce 
  public :: mpi_dist_cyclic_t

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

#ifdef WITH_MPI
    N = this%global_N / this%size
    ! Figure out if this node has any remaining elements
    if ( this%glob_2_rank(this%global_N) <= this%rank ) then
       N = N + ONE
    end if
#else
    ! Non MPI
    N = this%global_N
#endif
    
  end function local_N_

  function global_2_local_(this, global) result(local)
    class(mpi_dist_cyclic_t), intent(in) :: this
    integer(ii_), intent(in) :: global
    integer(ii_) :: local
    
#ifdef WITH_MPI
    ! Figure out if this node has this element
    if ( this%glob_2_rank(global) == this%rank ) then
       local = global / this%size
    else
       local = -ONE
    end if
#else
    ! Non MPI
    local = global
#endif
    
  end function global_2_local_

  function local_2_global_(this, local) result(global)
    class(mpi_dist_cyclic_t), intent(in) :: this
    integer(ii_), intent(in) :: local
    integer(ii_) :: global
    
#ifdef WITH_MPI
    global = this%size * local + this%rank
#else
    ! Non MPI
    global = local
#endif
    
  end function local_2_global_

  function global_2_rank_(this, global) result(rank)
    class(mpi_dist_cyclic_t), intent(in) :: this
    integer(ii_), intent(in) :: global
    integer(ii_) :: rank

#ifdef WITH_MPI
    rank = mod(global - ONE, this%global_N)
#else
    ! Non MPI
    rank = 0
#endif
    
  end function global_2_rank_

end module mpi_dist_cyclic
