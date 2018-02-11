module esl_sparse_matrix_m

  use esl_sparse_pattern_m, only: sparse_pattern_t
  use prec, only: dp

  implicit none

  public :: sparse_matrix_t

  !< A CSR sparse matrix in 1D
  !<
  !< A sparse matrix contains a reference to a sparse pattern.
  !< It is important when constructing this sparse matrix that the
  !< pattern is finalized before entry.
  !< If not, the program will stop.
  !<
  !< While this object is in use, the sparse pattern may *not* change
  !< elements. There is no way that the sparse matrix and/or sparse pattern
  !< has knowledge of the other one changing. Hence, we MUST keep them fixed.
  type sparse_matrix_t

    !< A pointer to the sparse pattern that defines this matrix
    type(sparse_pattern_t), pointer :: sp => null()
    real(dp), allocatable :: M(:)

  contains

    !< Initialize a sparse matrix
    procedure, public :: init => init_

    !< Print, to std-out information regarding this sparse pattern
    procedure, public :: print => print_

    !< Delete the sparse data (but not the sparse pattern)
    procedure, public :: delete => delete_

  end type sparse_matrix_t

contains

  subroutine init_(this, sparse_pattern)
    class(sparse_matrix_t), intent(inout) :: this
    class(sparse_pattern_t), intent(in), target :: sparse_pattern

    this%sp => sparse_pattern
    if ( .not. sparse_pattern%finalized ) then
      stop 'The sparse pattern has not been finalized before creating a sparse matrix'
    end if

    if ( allocated(this%M) ) then
      deallocate(this%M)
    end if
    allocate(this%M(this%sp%nt))

  end subroutine init_

  subroutine delete_(this, stat)
    class(sparse_matrix_t), intent(inout) :: this
    integer, intent(out), optional :: stat
    integer :: istat

    ! Ensure clean pointer
    nullify(this%sp)

    if ( present(stat) ) stat = 0
    if ( .not. allocated(this%M) ) return
    deallocate(this%M, stat=istat)
    if ( present(stat) ) stat = istat

  end subroutine delete_

  subroutine print_(this)
    class(sparse_matrix_t), intent(in) :: this
    write(*,'(a)') "<sparse_matrix>"
    call this%sp%print()
    write(*,'(a)') "</sparse_matrix>"
  end subroutine print_

end module esl_sparse_matrix_m





