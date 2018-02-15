!< The block cyclic distribution.
!<
!< The distribution is consecutive and each node gets
!< a block of elements in a cyclic manner.
!< Any remaining elements are added in a round-robin fashion in
!< block sizes.
!< This is equivalent to the scalapack distribution.
!< Note that, contrary to scalapack this forces the first
!< matrix processor to have rank 0.
module mpi_dist_block_cyclic_m

#ifdef WITH_MPI
  use mpi
#endif
  use mpi_dist_m

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
  public :: mpi_dist_block_cyclic_t

  integer(ii_), parameter :: ONE = 1_ii_

  type, extends(mpi_dist_t) :: mpi_dist_block_cyclic_t

    !< Block size of the cyclic distribution
    integer(ii_) :: block = ONE
    
  contains
    
    !< Initialize a new MPI distribution
    procedure, public :: init => init_
    
    !< Delete this object (it does NOT disconnect the communicator)
    procedure, public :: delete => delete_
    
    procedure, public :: glob_2_loc => global_2_local_
    procedure, public :: loc_2_glob => local_2_global_
    procedure, public :: glob_2_rank => global_2_rank_
    
  end type mpi_dist_block_cyclic_t

contains

  subroutine init_(this, comm, global_N, block)
    class(mpi_dist_block_cyclic_t), intent(inout) :: this
    !< MPI communicator
    integer(im_), intent(in) :: comm 
    !< Global number of elements
    integer(ii_), intent(in) :: global_N
    !< Block size of the distribution
    integer(ii_), intent(in) :: block

    call this%mpi_dist_t%init_(comm, global_N)
    this%block = block

  end subroutine init_

  subroutine delete_(this)
    class(mpi_dist_block_cyclic_t), intent(inout) :: this

    this%block = ONE
    call this%mpi_dist_t%delete_()

  end subroutine delete_

  function local_N_(this) result(N)
    class(mpi_dist_block_cyclic_t), intent(in) :: this
    integer(ii_) :: N

#ifdef WITH_MPI
    integer(ii_) :: iblock

    ! Abstracted from ScaLAPACK/numroc

    iblock = this%global_N / this%block
    ! For sure this node has this amount of elements
    N = (iblock / this%size) * this%block

    iblock = mod( iblock, this%size)
    if ( this%rank < iblock ) then
      ! less than the edeg
      N = N + this%block
    else if ( this%rank == iblock ) then
      ! edge rank (remaining elements
      N = N + mod( this%global_N, this%block)
    end if
#else
    ! Non MPI
    N = this%global_N
#endif
  end function local_N_

  function global_2_local_(this, global) result(local)
    class(mpi_dist_block_cyclic_t), intent(in) :: this
    integer(ii_), intent(in) :: global
    integer(ii_) :: local

#ifdef WITH_MPI
    ! Abstracted from ScaLAPACK/indxg2l
    local = this%block*( (global-ONE)/(this%block*this%size)) + mod(global-ONE,this%block)+ONE
#else
    ! Non MPI
    local = global
#endif

  end function global_2_local_

  function local_2_global_(this, local) result(global)
    class(mpi_dist_block_cyclic_t), intent(in) :: this
    integer(ii_), intent(in) :: local
    integer(ii_) :: global

#ifdef WITH_MPI
    ! Abstracted from ScaLAPACK/indxl2g
    global = this%size * this%block * ((local-ONE) / this%block) + &
        mod(local - ONE, this%block) + this%rank * this%block + ONE
#else
    ! Non MPI
    global = local
#endif

  end function local_2_global_

  function global_2_rank_(this, global) result(rank)
    class(mpi_dist_block_cyclic_t), intent(in) :: this
    integer(ii_), intent(in) :: global
    integer(ii_) :: rank

#ifdef WITH_MPI
    ! Abstracted from ScaLAPACK/indxg2p
    rank = mod( (global - ONE) / this%block, this%size)
#else
    ! Non MPI
    rank = 0
#endif

  end function global_2_rank_

end module mpi_dist_block_cyclic_m
