module sparse_pattern

  ! The default container integer kind
  use prec, only: ii_ => ip
  
  implicit none
  
  private

  ! Often used private variables
  integer(ii_), parameter :: ONE = 1_ii_
  integer(ii_), parameter :: ZERO = 0_ii_

  public :: sparse_pattern_t

  !< A CSR sparse pattern
  !<
  !< This is merely a container for a sparse pattern and can thus not
  !< be regarding as anything but rows and columns. There is no data
  !< stored in this type.
  !<
  !< An important feature for this sparse pattern is that it can be
  !< extended at any given time. I.e. one pre-allocates a "big" matrix
  !< with no sparse elements. The pool of elements can then be
  !< consecutively populated via the `this%add(r, c)` method.
  type sparse_pattern_t
     
     integer(ii_) :: nr = ZERO !< Number of rows in sparse pattern
     integer(ii_) :: nc = ZERO !< Number of columns in sparse pattern
     integer(ii_) :: nz = ZERO !< Number of non-zero elements in sparse pattern
     integer(ii_) :: nt = ZERO !< Maximum number of non-zero elements in sparse pattern
     
     logical :: finalized = .false. !< Parameter to check whether the sparse pattern is finalized (nz == nt)

     !< 1-based pointers for each row.
     !<   column(rptr(2))
     !< is the first column in the sparse pattern for row 2.
     !<   column(rptr(2) + nrow(2) - 1)
     !< is the last column in the sparse pattern for row 2.
     integer(ii_), allocatable :: rptr(:)
     
     !< Number of entries per row (size() == nr), sum(nrow) == nz
     !< This is a necessity only when constructing the matrix.
     !< It however, has some added benefits when one wishes to figure
     !< out the edges of a given node.
     integer(ii_), allocatable :: nrow(:)
     
     !< Column indices.
     integer(ii_), allocatable :: column(:)
     
   contains

     procedure, public :: init => init_dim_

     !< True if the sparse pattern is initialized
     procedure, public :: initialized => initialized_

     !< Query number of nonzero elements currently in the sparse pattern
     procedure, public :: nonzeros => nonzeros_
     !< Query maximum number of nonzero elements that this sparse pattern may contain
     procedure, public :: max_nonzeros => max_nonzeros_

     !< Number of rows in the sparse pattern 
     procedure, public :: rows => rows_
     !< Number of columns in the sparse pattern 
     procedure, public :: columns => columns_
     
     !< Directly query the sparse index for a given `row, column`
     procedure, public :: index => index_

     !< Add a new row and column to the sparse pattern
     procedure, public :: add => add_

     !< Finalize the sparse pattern such that the pointers and elements
     procedure, public :: finalize => finalize_

     !< Print, to std-out information regarding this sparse pattern
     procedure, public :: print => print_

     !< Delete this sparse object
     procedure, public :: delete => delete_
     
  end type sparse_pattern_t
  
contains

  function initialized_(this) result(initd)
    class(sparse_pattern_t), intent(in) :: this
    logical :: initd
    initd = allocated(this%column)
  end function initialized_

  function nonzeros_(this) result(nnzs)
    class(sparse_pattern_t), intent(in) :: this
    integer(ii_) :: nnzs
    nnzs = this%nz
  end function nonzeros_

  function max_nonzeros_(this) result(nt)
    class(sparse_pattern_t), intent(in) :: this
    integer(ii_) :: nt
    nt = this%nt
  end function max_nonzeros_
  
  function rows_(this) result(rows)
    class(sparse_pattern_t), intent(in) :: this
    integer(ii_) :: rows
    rows = this%nr
  end function rows_
  
  function columns_(this) result(cols)
    class(sparse_pattern_t), intent(in) :: this
    integer(ii_) :: cols
    cols = this%nc
  end function columns_
  
  subroutine print_(this)
    class(sparse_pattern_t), intent(in) :: this
    write(*,'(3(a,i0),a)') "<sparse_pattern rows=", this%nr, &
         ", columns=", this%nc, ", non-zeros=",this%nz, "/>"
  end subroutine print_
  
  subroutine delete_(this, stat)
    class(sparse_pattern_t), intent(inout) :: this
    integer, intent(out), optional :: stat
    integer :: istat, max_stat
    this%nr = 0
    this%nc = 0
    this%nz = 0
    this%nt = 0
    this%finalized = .false.
    
    if ( present(stat) ) stat = 0
    if ( .not. allocated(this%column) ) return
    
    deallocate(this%rptr, stat=istat)
    max_stat = istat
    deallocate(this%column, stat=stat)
    if ( istat /= 0 ) max_stat = istat
    deallocate(this%nrow, stat=stat)
    if ( istat /= 0 ) max_stat = istat

    if ( present(stat) ) stat = max_stat
    
  end subroutine delete_
  
  subroutine init_dim_(this, nr, nc, nt, np)
    class(sparse_pattern_t), intent(inout) :: this
    integer(ii_), intent(in) :: nr, nc
    integer(ii_), intent(in), optional :: nt, np
    integer(ii_) :: i, lnp

    ! store the pattern size
    this%nr = nr
    this%nc = nc
    ! store the initial size of the sparse pattern
    if ( present(nt) ) then
       this%nt = nt
    else if ( present(np) ) then
       this%nt = np * nr
    else
       ! np is defaulting to 10
       this%nt = nr * 10_ii_
    end if

    ! initial number of non-zero elements
    this%nz = ZERO
    ! this object is *not* finalized
    this%finalized = .false.

    ! Allocate data-structure
    allocate(this%rptr(nr+ONE))
    allocate(this%nrow(nr))
    this%nrow = ZERO
    allocate(this%column(this%nt))

    ! pre-define pointers to sparse pattern
    lnp = this%nt / this%nr
    this%rptr(ONE) = ONE
    do i = ONE , nr
       this%rptr(i) = this%rptr(i-ONE) + np
    end do
    this%rptr(nr+ONE) = this%rptr(nr) + ONE
    
  end subroutine init_dim_

  !< Figure out index of row, column for the sparse pattern. If the element is zero (non-existing) we return -1
  pure function index_(this, ir, ic) result(idx)
    class(sparse_pattern_t), intent(in) :: this
    integer(ii_), intent(in) :: ir
    integer(ii_), intent(in) :: ic
    integer(ii_) :: idx
    do idx = this%rptr(ir) , this%rptr(ir) + this%nrow(ir) - ONE
       if ( this%column(idx) == ic ) return
    end do
    idx = -ONE
  end function index_

  subroutine finalize_(this)
    class(sparse_pattern_t), intent(inout) :: this
    integer(ii_) :: ind, ir, rind, first_ind
    integer(ii_), allocatable :: col(:)
    
    ! Current "new" index point
    ind = ONE
    
    ! Loop rows
    do ir = ONE, this%nr
       
       ! Loop over sparse elements
       first_ind = ind
       do rind = this%rptr(ir), this%rptr(ir) + this%nrow(ir) - ONE
          this%column(ind) = this%column(rind)
          ind = ind + ONE
       end do
       
       this%rptr(ir) = first_ind
       
    end do

    allocate(col(ind - ONE))
    col(:) = this%column(:ind - ONE)
    deallocate(this%column)
    allocate(this%column(ind-ONE))
    this%column = col
    this%nt = ind - ONE

    ! Update last pointer
    this%rptr(this%nr) = ind
    
    this%finalized = .true.
    
  end subroutine finalize_
  
  subroutine add_(this, ir, ic, ind)
    class(sparse_pattern_t), intent(inout) :: this
    integer(ii_), intent(in) :: ir, ic
    integer(ii_), intent(out), optional :: ind
    ! Used in sub-routines
    integer(ii_) :: lind

    ! first try simple add
    if ( simple_add() ) then
       if ( present(ind) ) ind = lind
       return
    end if

    ! We have to do the complex thing... :(
    ! This means heavy work due to deallocation, allocation, etc.
    call complex_add()
    if ( present(ind) ) ind = lind
    
  contains

    function simple_add() result(add)
      logical :: add

      ! rptr(ir) + nrow(ir) == next empty place in sparse pattern
      ! So if the following pointer is placed further, then proceed
      lind = this%rptr(ir) + this%nrow(ir)
      ! Check whether there is room
      add = this%rptr(ir + ONE) - lind > ZERO

      if ( .not. add ) return

      this%column(lind) = ic
      this%nrow(ir) = this%nrow(ir) + ONE

    end function simple_add

    subroutine complex_add()
      ! A copy of cols before we insert it.
      integer(ii_), allocatable :: cols(:)

      stop 'currently not implemented'
      
    end subroutine complex_add

  end subroutine add_
  
end module sparse_pattern
