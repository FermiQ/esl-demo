!< This module implements a generic MPI based distribution
!<
!< This type should be an inherited class to ensure that
!< we can implement different distributions.
!<
!< The distribution function may be used to figure out
!< localities of linear quantities, i.e. 1D problems.
module mpi_dist

#ifdef WITH_MPI
  use mpi
#endif

  ! MPI relies on standard integer sizes
  ! However, if one wishes to build this application
  ! with default long integers, we have to be careful.
  use prec, only: im_ => ip

  ! To increase precision for the container
  ! simply change the pointed to variable.
  use prec, only: ii_ => ip

  implicit none

  private

  public :: mpi_dist_t

  ! There are problems with having abstract interfaces and
  ! methods callable from sub-classed types.
  ! Hence we must simply force the implementor of new distributions
  ! to implement:
  !
  !   glob_2_loc
  !   loc_2_glob
  !   glob_2_rank
  !   local_N
  type :: mpi_dist_t

     !< Associated communicator with this distribution.
     integer(im_) :: comm = -1

     ! Note that rank and size are kept in the same precision
     ! as global_N because it allows "faster" arithmetic due to
     ! less casting. Or at least ideally...

     !< The MPI rank (MPI_Comm_rank(rank)) (i.e. 0-based)
     integer(ii_) :: rank = 0

     !< The MPI communicator size (MPI_Comm_size(size))
     integer(ii_) :: size = 1

     !< Total number of elements to distribute
     integer(ii_) :: global_N = 0

   contains

     ! Sadly creating public procedures (or generic procedures)
     ! forces the interfaces to be *exactly* the same.
     ! Perhaps I am missing something?
     ! I.e. sub-classed types should call %mpi_dist_t%init_(..)
     ! for correct interaction.
     ! This may propagate further by introducing problems when
     ! inheriting sub-classed distributions where the interface
     ! may be required to change.

     !< Initialize a new MPI distribution
     procedure, public, non_overridable :: init_

     !< Delete this object (it does NOT disconnect the communicator)
     procedure, public, non_overridable :: delete_

  end type mpi_dist_t

contains

  subroutine init_(this, comm, global_N)
    class(mpi_dist_t), intent(inout) :: this
    !< MPI communicator
    integer(im_), intent(in) :: comm 
    !< Global number of elements
    integer(ii_), intent(in) :: global_N

    integer(im_) :: err

    ! Figure out the MPI things
    this%comm = comm
#ifdef WITH_MPI
    call MPI_Comm_rank(comm, this%rank, err)
    call MPI_Comm_rank(comm, this%size, err)
#else
    this%rank = 0
    this%size = 1
#endif
    this%global_N = global_N

  end subroutine init_

  subroutine delete_(this)
    class(mpi_dist_t), intent(inout) :: this

    this%comm = -1_ii_
    this%rank = 0_ii_
    this%size = 1_ii_
    this%global_N = 0_ii_

  end subroutine delete_

end module mpi_dist
