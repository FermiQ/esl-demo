!< Stub module to allow ELSI calls without having elsi
!<
!< This also implements the smearing methods
module elsi
  use prec, only: ip, dp

  implicit none

  ! Everything is public (except precisions)
  public
  private :: ip, dp

  type elsi_handle
    real(dp) :: nel
    integer(ip) :: method
    real(dp) :: width
  end type elsi_handle

contains

  subroutine elsi_init(eh, solver, mode, mat_format, nbasis, nel, nstate)
    type(elsi_handle), intent(inout) :: eh
    real(dp), intent(in) :: nel
    integer(ip), intent(in) :: solver, mode, mat_format, nbasis, nstate
    eh%nel = nel
  end subroutine elsi_init

  subroutine elsi_finalize(eh)
    type(elsi_handle), intent(inout) :: eh
  end subroutine elsi_finalize

  subroutine elsi_set_csc(eh, global_nnz, local_nnz, nrow, col_idx, row_ptr)
    type(elsi_handle), intent(inout) :: eh
    integer(ip), intent(in) :: global_nnz, local_nnz, nrow
    integer(ip), intent(in) :: col_idx(:), row_ptr(:)
  end subroutine elsi_set_csc

  subroutine elsi_dm_real_sparse(eh, h, s, d, energy)
    type(elsi_handle), intent(in) :: eh
    real(dp), intent(inout) :: h(:), s(:), d(:)
    real(dp), intent(out) :: energy
  end subroutine elsi_dm_real_sparse

  subroutine elsi_get_edm_real_sparse(eh, edm)
    type(elsi_handle), intent(in) :: eh
    real(dp), intent(inout) :: edm(:)
  end subroutine elsi_get_edm_real_sparse

  subroutine elsi_compute_mu_and_occ(eh, nel, nstate, nspin, nkpt, evals, &
      occ, k_weights, mu)
    type(elsi_handle), intent(in) :: eh
    real(dp), intent(in) :: nel
    integer(ip), intent(in) :: nstate, nspin, nkpt
    real(dp), dimension(nstate,nspin,nkpt), intent(in) :: evals
    real(dp), dimension(nstate,nspin,nkpt), intent(out) :: occ
    real(dp), dimension(nkpt), intent(in) :: k_weights
    real(dp), intent(out) :: mu

    real(dp) :: e_min, e_max, spin_degen
    real(dp) :: mu_min, mu_max
    real(dp) :: diff_min, diff_max, diff

    e_min = minval(evals)
    e_max = maxval(evals)

    ! Initialize
    mu_min = e_min
    mu_max = e_max
    if ( mu_max - mu_min < 1.e-6_dp ) then
      mu_max = mu_max + 1._dp
    end if

    ! Initialize
    occ(:,:,:) = 0._dp

    call next_step()
    !print *, 'occ', occ

    ! Currently we have no fall-back
    do while ( diff_min * diff_max > 0._dp )

      mu_min = mu_min - 0.5_dp * abs(e_max - e_min)
      mu_max = mu_max + 0.5_dp * abs(e_max - e_min)

      call next_step()

    end do

    ! Now calculate the Fermi-level
    if ( abs(diff_min) < 1.e-6_dp ) then
      mu = mu_min
      call single_step(mu, diff)
      return
    else if ( abs(diff_max) < 1.e-6_dp ) then
      mu = mu_max
      call single_step(mu, diff)
      return
    end if

    ! TODO clean up to really have a correctly behaving bisection.

    do
      mu = (mu_min + mu_max) * 0.5_dp
      call single_step(mu, diff)
      if ( abs(diff) < 1.e-6_dp ) then
        return
      else if ( diff < 0._dp ) then
        mu_min = mu
      else ! diff > 0._dp
        mu_max = mu
      end if
    end do

  contains

    subroutine next_step()

      select case ( eh%method )
      case ( 1 ) ! Fermi-Dirac
        spin_degen = 3 - nspin
        call fermi_dirac_cycle(mu_min, diff_min)
        call fermi_dirac_cycle(mu_max, diff_max)
      case default
        print *, 'ELSI stub does not implement others than Fermi-Dirac distributions'
      end select

    end subroutine next_step

    subroutine single_step(mu, diff)
      real(dp), intent(in) :: mu
      real(dp), intent(out) :: diff

      select case ( eh%method )
      case ( 1 ) ! Fermi-Dirac
        spin_degen = 3 - nspin
        call fermi_dirac_cycle(mu, diff)
      case default
        print *, 'ELSI stub does not implement others than Fermi-Dirac distributions'
      end select

    end subroutine single_step

    subroutine fermi_dirac_cycle(mu, diff)
      real(dp), intent(in) :: mu
      real(dp), intent(out) :: diff


      integer :: ikpt, ispin, istate

      real(dp) :: max_exp, invert_width, arg

      invert_width = 1._dp/eh%width
      diff = 0._dp
      max_exp = maxexponent(mu)*log(2.0_dp)

      do ikpt = 1, nkpt
        do ispin = 1, nspin
          do istate = 1, nstate
            arg = (evals(istate,ispin,ikpt) - mu) * invert_width
            if ( arg < max_exp ) then
              occ(istate,ispin,ikpt) = spin_degen / (1._dp + exp(arg))
              diff = diff + occ(istate,ispin,ikpt) * k_weights(ikpt)
            else
              occ(istate,ispin,ikpt) = 0._dp
            end if
          end do
        end do
      end do

      diff = diff - eh%nel

    end subroutine fermi_dirac_cycle
    
  end subroutine elsi_compute_mu_and_occ

  subroutine elsi_set_mu_broaden_scheme(eh, method)
    type(elsi_handle), intent(inout) :: eh
    integer(ip), intent(in) :: method
    eh%method = method
  end subroutine elsi_set_mu_broaden_scheme

  subroutine elsi_set_mu_broaden_width(eh, width)
    type(elsi_handle), intent(inout) :: eh
    real(dp), intent(in) :: width
    eh%width = width
  end subroutine elsi_set_mu_broaden_width

end module elsi
