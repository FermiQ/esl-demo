module esl_system_m

  use prec, only : dp,ip
  use esl_numeric_m, only : matr3inv

  use fdf, only : block_fdf, fdf_get, fdf_block,fdf_defined, &
       parsed_line, fdf_breals, fdf_bline, fdf_bnames

  ! Sparse pattern for LO
  use esl_sparse_pattern_m, only: sparse_pattern_t
  ! Sparse matrix for LO (only real(dp))
  use esl_sparse_matrix_m, only: sparse_matrix_t

  use esl_basis_m
  use esl_smear_m
  use esl_states_m
  use esl_species_m
  use esl_geometry_m
  
  implicit none
  private

  public ::                &
       system_t

  !Data structure for the system
  type system_t
    type(basis_t)    :: basis
    type(geometry_t) :: geo
    
    ! TODO decide whether the system should be inherited for the LO/PW
    ! case. It may make the system a lot easier to figure out.
    ! However, it will prohibit switching from PW/LO -> LO/PW within the same
    ! calculation.

    ! LO dependent variables
    type(sparse_pattern_t):: sparse_pattern
    type(sparse_matrix_t) :: S ! always 1D
    type(sparse_matrix_t), allocatable :: H(:) ! one per spin
    type(sparse_matrix_t), allocatable :: DM(:) ! one per spin

    real(dp) :: nElectrons
  contains
    private
    procedure, public :: init
    procedure, public :: summary
    final  :: cleanup
  end type system_t

contains

  !Initialize the physical system
  !----------------------------------------------------
  subroutine init(sys)
    class(system_t) :: sys

    integer :: nstates, nspin

    call sys%geo%init()
    
    call sys%basis%init(sys%geo)

    call sys%summary()

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(sys)
    type(system_t) :: sys

  end subroutine cleanup

  !Summary
  !----------------------------------------------------
  subroutine summary(sys)
    use yaml_output
    class(system_t) :: sys

    call yaml_mapping_open("System")
    call sys%geo%summary()
    call sys%basis%summary()
    call yaml_mapping_close()

  end subroutine summary

end module esl_system_m
