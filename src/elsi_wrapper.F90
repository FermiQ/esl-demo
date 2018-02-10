module elsi_wrapper_esl
  use prec
  use elsi
  implicit none
  private

  public :: elsi_t
  public :: elsi_calc_dm            !< Compute density matrix
  !public :: elsi_calc_ev            !< Compute eigensolutions
  public :: elsi_calc_fermi_and_occ !< Compute Fermi level and occupations

  !Data structure for ELSI
  !TODO: Wait for the sparse matrix type
  type elsi_t
     type(elsi_handle)   :: e_h
     type(sparse_matrix) :: H
     type(sparse_matrix) :: S
     type(sparse_matrix) :: D
     real(kind=dp)       :: energy
   contains
     private
     procedure, public :: init
     final :: cleanup
  end type elsi_t

  integer(kind=ip), parameter :: ELPA        = 1
  integer(kind=ip), parameter :: LIBOMM      = 2
  integer(kind=ip), parameter :: PEXSI       = 3
  integer(kind=ip), parameter :: SIPS        = 5
  integer(kind=ip), parameter :: SINGLE_PROC = 0
  integer(kind=ip), parameter :: MULTI_PROC  = 1
  integer(kind=ip), parameter :: PEXSI_CSC   = 1

contains

  !Initialize ELSI
  !----------------------------------------------------
  subroutine init(this)
    class(elsi_t) :: this

    !Initialize an ELSI handle
    !TODO: Feed in system info
    !call elsi_init(this%e_h,ELPA,MULTI_PROC,PEXSI_CSC,n_basis,n_electron,n_state)

    !Initialize ELSI MPI
    !TODO: Feed in whatever communicator
    !call elsi_set_mpi(this%e_h,mpi_comm)

    !Initialize ELSI SPARSITY PATTERN
    !TODO: Feed in sparse matrix info
    !call elsi_set_csc(this%e_h,global_nnz,local_nnz,local_ncol,row_idx,col_ptr)

  end subroutine init

  !Finalize ELSI
  !----------------------------------------------------
  subroutine cleanup(this)
    type(elsi_t) :: this

    call elsi_finalize(this%e_h)

  end subroutine cleanup

  !Compute density matrix by ELSI
  !----------------------------------------------------
  subroutine elsi_calc_dm(this)
    type(elsi_t), intent(inout) :: this

    !Only real value for now
    call elsi_dm_real(this%e_h,this%H%dat,this%S%dat,this%D%dat,this%energy)

  end subroutine elsi_calc_dm

  !Compute density matrix by ELSI
  !----------------------------------------------------
  subroutine elsi_calc_fermi_and_occ(this, n_electron, n_state, n_spin, n_kpt, &
       & evals, occ_nums, k_weights, fermi)
    type(elsi_t), intent(inout) :: this
    integer, intent(in) :: n_electron, n_state, n_spin, n_kpt
    real(dp), dimension(n_state,n_spin,n_kpt), intent(in) :: evals
    real(dp), dimension(n_state,n_spin,n_kpt), intent(in) :: occ_nums
    real(dp), dimension(n_kpt), intent(in) :: k_weights
    real(dp), intent(out) :: fermi

    !TODO: Wait for the input info
    call elsi_compute_mu_and_occ(this%e_h,n_electron,n_state,n_spin,n_kpt,&         
         k_weights,evals,occ_nums,fermi)

  end subroutine elsi_calc_fermi_and_occ

end module elsi_wrapper_esl
