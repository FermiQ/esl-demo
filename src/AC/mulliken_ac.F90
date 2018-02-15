!< Calculation of the Mulliken charges and printing them out to YAML
module esl_mulliken_ac_m

  use prec
  use esl_basis_ac_m
  use esl_sparse_pattern_m
  use esl_sparse_matrix_m

  implicit none

  private

  public :: mulliken_ac_summary

contains

  !< Calculate and print-out the Mulliken charges
  subroutine mulliken_ac_summary(basis, S, DM)
    use yaml_output

    type(basis_ac_t), intent(in) :: basis
    type(sparse_matrix_t), intent(in) :: S, DM

    ! Local variables for calculating the Mulliken charges
    integer :: ia, io, ind
    integer :: io1, io2

    type(sparse_pattern_t), pointer :: sp
    
    ! Actual Mulliken charges per site
    real(dp), allocatable :: M(:), F(:)
    character(len=10) :: str

    ! Retrieve sparse pattern
    sp => S%sp
    
    allocate(M(basis%n_site))
    allocate(F(basis%n_orbital))

    ! Calculate the total number charge, per site
    do ia = 1 , basis%n_site
      
      M(ia) = 0._dp
      do io = basis%site_orbital_start(ia), basis%site_orbital_start(ia + 1) - 1

        F(io) = 0._dp
        do ind = sp%rptr(io), sp%rptr(io) + sp%nrow(io) - 1

          F(io) = F(io) + S%M(ind) * DM%M(ind)
          
        end do

        M(ia) = M(ia) + F(io)

      end do

    end do

    ! Produce YAML output
    call yaml_mapping_open("Mulliken")
    call yaml_map('Total', sum(M))
    call yaml_comment("Q", hfill = "-")
    do ia = 1, basis%n_site

      write(str, '(i0)') ia
      call yaml_mapping_open(trim(str))
      !str = 'StillNotCreated' !basis%species(basis%site_state_idx(ia))%label
      !call yaml_map('Label', trim(str))
      call yaml_map('Sum', M(ia))
      io1 = basis%site_orbital_start(ia)
      io2 = basis%site_orbital_start(ia + 1) - 1
      call yaml_map('Individual', F(io1:io2))
      call yaml_mapping_close()
      
    end do

    call yaml_mapping_close()

    deallocate(M, F)
    
  end subroutine mulliken_ac_summary

end module esl_mulliken_ac_m

