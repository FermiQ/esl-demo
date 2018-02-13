! Copyright (c) 2015-2017, the ELSI team. All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
!  * Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
!  * Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
!  * Neither the name of the "ELectronic Structure Infrastructure" project nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!>
!! This module contains subroutines to solve an eigenproblem using a reverse
!! communication interface.
!!
module ELSI_RCI_OMM

    use ELSI_RCI_BLAS
    use ELSI_RCI_CONSTANTS
    use ELSI_RCI_DATATYPE
    use ELSI_RCI_PRECISION, only:i4

    implicit none

    public :: rci_omm

    ! Matrix ID of size m by n
    !&<
    integer(kind=i4), parameter :: MID_C    = 1
    integer(kind=i4), parameter :: MID_HC   = 2
    integer(kind=i4), parameter :: MID_SC   = 3
    integer(kind=i4), parameter :: MID_G    = 4
    integer(kind=i4), parameter :: MID_PG   = 5
    integer(kind=i4), parameter :: MID_HG   = 6
    integer(kind=i4), parameter :: MID_SG   = 7
    integer(kind=i4), parameter :: MID_D    = 8
    integer(kind=i4), parameter :: MID_Gp   = 9
    integer(kind=i4), parameter :: MID_PGp  = 10
    integer(kind=i4), parameter :: MID_WORK = 11
    !&>

    ! Matrix ID of size n by n
    !&<
    integer(kind=i4), parameter :: MID_HW   = 21
    integer(kind=i4), parameter :: MID_SW   = 22
    integer(kind=i4), parameter :: MID_HWd  = 23
    integer(kind=i4), parameter :: MID_SWd  = 24
    integer(kind=i4), parameter :: MID_HWdd = 25
    integer(kind=i4), parameter :: MID_SWdd = 26
    integer(kind=i4), parameter :: MID_HWdT = 27
    integer(kind=i4), parameter :: MID_SWdT = 28
    !&>

    ! Stage ID
    !&<
    integer(kind=i4), parameter :: SID_INIT    = 0
    integer(kind=i4), parameter :: SID_COEFF   = 100
    integer(kind=i4), parameter :: SID_ITER    = 200
    integer(kind=i4), parameter :: SID_UPDATE  = 300
    integer(kind=i4), parameter :: SID_LS_FAIL = 400
    integer(kind=i4), parameter :: SID_LS_CONV = 500
    !&>

contains

    subroutine rci_omm(ijob, iS, task, result_in, m, n, e_min, &
                       cg_tol, max_iter, long_out)

        implicit none

        !**** INPUT ***********************************!

        logical, intent(in) :: long_out ! print detailed output?

        integer, intent(in) :: m ! size of basis
        integer, intent(in) :: n ! number of occupied states

        integer, intent(in)  :: max_iter
        real(r8), intent(in) :: cg_tol ! convergence tolerance of CG
        real(r8), intent(in) :: result_in(:)
        ! minimization (if negative, default of 1.0d-9 is used)

        !**** OUTPUT **********************************!
        ! OMM functional energy (spin degeneracy *not* included)
        real(r8), intent(out) :: e_min
        integer, intent(out) :: task

        !**** INOUT ***********************************!
        integer, intent(inout) :: ijob ! job id
        type(rci_instr), intent(inout) :: iS

        !**** LOCAL ***********************************!
        logical, save :: ls_conv ! line search converged?
        logical, save :: ls_fail ! line search failed?

        integer, save :: icg ! CG step num.
        integer, save :: mpi_size, mpi_rank

        real(r8), save :: cg_tol_internal
        real(r8), save :: lambda ! CG step size
        real(r8), save :: lambda_d ! lambda denominator
        real(r8), save :: lambda_n ! lambda numerator
        real(r8), save :: e_min_save
        real(r8), save :: e_min_old ! OMM func energy at previous CG step
        real(r8), save :: coeff(0:4) ! coeffs. of the quartic equation
        real(r8), save :: x_min ! line search position of minimum
        real(r8), save :: e_diff

        real(r8), save :: TrH
        real(r8), save :: TrHS
        real(r8), save :: TrHd
        real(r8), save :: TrHdS
        real(r8), save :: TrHSd
        real(r8), save :: TrHdd
        real(r8), save :: TrHddS
        real(r8), save :: TrHSdd
        real(r8), save :: TrHdSd
        real(r8), save :: TrHdSdT
        real(r8), save :: TrHddSd
        real(r8), save :: TrHdSdd
        real(r8), save :: TrHddSdd

        logical, save :: conv
        !**********************************************!

        ! if this is the first SCF step, then we need to initialize the WF
        ! coeffs. matrix with random numbers between -0.5 and 0.5 (normalize
        ! at the end to avoid instabilities), unless we are reading them
        ! from file

        ! first we calculate the energy and gradient for our initial guess,
        ! with the following steps:
        ! -calculate the hamiltonian in WF basis: HW=C^T*H*C
        ! -- HC = H*C
        if (ijob <= SID_INIT) then
            ijob = SID_INIT + 1
            if (cg_tol < 0.0_r8) then
                cg_tol_internal = 1.0d-9
            else
                cg_tol_internal = cg_tol
            end if
            icg = 0
            lambda = 0.0_r8
            call rci_h_multi(iS, task, 'N', m, n, MID_C, MID_HC)
            return
        end if
        ! -- HW = C'*HC
        if (ijob == SID_INIT + 1) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_C, MID_HC, MID_HW)
            ijob = ijob + 1
            return
        end if

        ! -calculate the overlap matrix in WF basis: SW=C^T*S*C
        ! -- SC = S*C
        if (ijob == SID_INIT + 2) then
            call rci_s_multi(iS, task, 'N', m, n, MID_C, MID_SC)
            ijob = ijob + 1
            return
        end if
        ! -- SW = C'*SC
        if (ijob == SID_INIT + 3) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_C, MID_SC, MID_SW)
            ijob = ijob + 1
            return
        end if

        ! -calculate the gradient: G=2*(2*HC-SC*HW-HC*SW)
        ! Note: It is required that G = 0 as the initial
        ! -- G = 4*HC
        if (ijob == SID_INIT + 4) then
            call rci_axpy(iS, task, m, n, 4.0_r8, MID_HC, MID_G)
            ijob = ijob + 1
            return
        end if
        ! -- G = -2 HC*SW + G
        if (ijob == SID_INIT + 5) then
            call rci_gemm(iS, task, 'N', 'N', m, n, n, -2.0_r8, 1.0_r8, &
                          MID_HC, MID_SW, MID_G)
            ijob = ijob + 1
            return
        end if
        ! -- G = -2 SC*HW + G
        if (ijob == SID_INIT + 6) then
            call rci_gemm(iS, task, 'N', 'N', m, n, n, -2.0_r8, 1.0_r8, &
                          MID_SC, MID_HW, MID_G)
            ijob = ijob + 1
            return
        end if

        ! -calculate the preconditioned gradient by premultiplying G by P
        ! -- PG = P*G
        if (ijob == SID_INIT + 7) then
            call rci_p_multi(iS, task, 'N', m, n, MID_G, MID_PG)
            ijob = ijob + 1
            return
        end if

        ! - calculate HWd = G'*H*C = PG'*HC
        if (ijob == SID_INIT + 8) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_PG, MID_HC, MID_HWd)
            ijob = ijob + 1
            return
        end if

        ! - calculate SWd = G'*S*C = PG'*SC
        if (ijob == SID_INIT + 9) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_PG, MID_SC, MID_SWd)
            ijob = ijob + 1
            return
        end if

        ! - calculate HWdd = G'*H*G = PG'*(H*PG)
        ! -- HG = H*PG
        if (ijob == SID_INIT + 10) then
            call rci_h_multi(iS, task, 'N', m, n, MID_PG, MID_HG)
            ijob = ijob + 1
            return
        end if
        ! -- HWdd = PG'*HG
        if (ijob == SID_INIT + 11) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_PG, MID_HG, MID_HWdd)
            ijob = ijob + 1
            return
        end if

        !  SWdd=G^T*S*G
        ! - calculate SWdd = G'*S*G = PG'*(S*PG)
        ! -- SG = S*PG
        if (ijob == SID_INIT + 12) then
            call rci_s_multi(iS, task, 'N', m, n, MID_PG, MID_SG)
            ijob = ijob + 1
            return
        end if
        ! -- SWdd = PG'*SG
        if (ijob == SID_INIT + 13) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_PG, MID_SG, MID_SWdd)
            ijob = SID_COEFF
            return
        end if

        ! calculate coefficients
        ! - TrH = trace(HW)
        if (ijob == SID_COEFF) then
            call rci_trace(iS, task, n, n, MID_HW)
            ijob = ijob + 1
            return
        end if
        ! - TrHS = trace(HW*SW)
        if (ijob == SID_COEFF + 1) then
            TrH = result_in(1)
            call rci_dot(iS, task, n, n, MID_HW, MID_SW)
            ijob = ijob + 1
            return
        end if
        ! - TrHd = trace(HWd)
        if (ijob == SID_COEFF + 2) then
            TrHS = result_in(1)
            call rci_trace(iS, task, n, n, MID_HWd)
            ijob = ijob + 1
            return
        end if
        ! - TrHdS = trace(HWd*SW)
        if (ijob == SID_COEFF + 3) then
            TrHd = result_in(1)
            call rci_dot(iS, task, n, n, MID_HWd, MID_SW)
            ijob = ijob + 1
            return
        end if
        ! - TrHSd = trace(HW*SWd)
        if (ijob == SID_COEFF + 4) then
            TrHdS = result_in(1)
            call rci_dot(iS, task, n, n, MID_HW, MID_SWd)
            ijob = ijob + 1
            return
        end if
        ! - TrHdd = trace(HWdd)
        if (ijob == SID_COEFF + 5) then
            TrHSd = result_in(1)
            call rci_trace(iS, task, n, n, MID_HWdd)
            ijob = ijob + 1
            return
        end if
        ! - TrHddS = trace(HWdd*SW)
        if (ijob == SID_COEFF + 6) then
            TrHdd = result_in(1)
            call rci_dot(iS, task, n, n, MID_HWdd, MID_SW)
            ijob = ijob + 1
            return
        end if
        ! - TrHSdd = trace(HW*SWdd)
        if (ijob == SID_COEFF + 7) then
            TrHddS = result_in(1)
            call rci_dot(iS, task, n, n, MID_HW, MID_SWdd)
            ijob = ijob + 1
            return
        end if
        ! - TrHdSdT = trace(HWd*SWdT)
        ! - SWdT = SWd'
        if (ijob == SID_COEFF + 8) then
            TrHSdd = result_in(1)
            call rci_copy(iS, task, 'C', n, n, MID_SWd, MID_SWdT)
            ijob = ijob + 1
            return
        end if
        if (ijob == SID_COEFF + 9) then
            call rci_dot(iS, task, n, n, MID_HWd, MID_SWdT)
            ijob = ijob + 1
            return
        end if
        ! - TrHdSd = trace(HWd*SWd)
        if (ijob == SID_COEFF + 10) then
            TrHdSdT = result_in(1)
            call rci_dot(iS, task, n, n, MID_HWd, MID_SWd)
            ijob = ijob + 1
            return
        end if
        ! - TrHddSd = trace(HWdd*SWd)
        if (ijob == SID_COEFF + 11) then
            TrHdSd = result_in(1)
            call rci_dot(iS, task, n, n, MID_HWdd, MID_SWd)
            ijob = ijob + 1
            return
        end if
        ! - TrHdSdd = trace(HWd*SWdd)
        if (ijob == SID_COEFF + 12) then
            TrHddSd = result_in(1)
            call rci_dot(iS, task, n, n, MID_HWd, MID_SWdd)
            ijob = ijob + 1
            return
        end if
        ! - TrHddSdd = trace(HWdd*SWdd)
        if (ijob == SID_COEFF + 13) then
            TrHdSdd = result_in(1)
            call rci_dot(iS, task, n, n, MID_HWdd, MID_SWdd)
            ijob = ijob + 1
            return
        end if
        ! - calculate coeff
        if (ijob == SID_COEFF + 14) then
            call rci_null(iS, task)
            TrHddSdd = result_in(1)
            coeff(0) = 2.0_r8*TrH - TrHS
            coeff(1) = 2.0_r8*(2.0_r8*TrHd - TrHdS - TrHSd)
            coeff(2) = 2.0_r8*(TrHdd - TrHdSd - TrHdSdT) - TrHddS - TrHSdd
            coeff(3) = -2.0_r8*(TrHddSd + TrHdSdd)
            coeff(4) = -TrHddSdd
            e_min = coeff(0)
            e_min_save = coeff(0)
            if (icg > 0) then
                ijob = SID_UPDATE
            else
                ijob = SID_ITER
            end if
            return
        end if

        ! this is the main loop of the CG algorithm. We perform a series
        ! of line minimizations, with the gradient G at each new step being
        ! modified to obtain the search direction D

        ! - D = PG
        if (ijob == SID_ITER) then
            call rci_copy(iS, task, 'N', m, n, MID_PG, MID_D)
            ijob = ijob + 1
            return
        end if
        ! - G_p = G
        if (ijob == SID_ITER + 1) then
            call rci_copy(iS, task, 'N', m, n, MID_G, MID_Gp)
            ijob = ijob + 1
            return
        end if
        ! - PG_p = PG
        if (ijob == SID_ITER + 2) then
            call rci_copy(iS, task, 'N', m, n, MID_PG, MID_PGp)
            ijob = ijob + 1
            return
        end if
        ! - e_min_old = e_min
        if (ijob == SID_ITER + 3) then
            call rci_null(iS, task)
            e_min_old = e_min_save
            if (icg <= 0) then
                ijob = SID_UPDATE
            else
                ijob = ijob + 1
            end if
            return
        end if

        ! if this is not the first CG step, we have to recalculate HWd,
        ! SWd, HWdd, SWdd, and the coeffs.
        ! if icg > 0
        ! HWd = D'*HC
        if (ijob <= SID_ITER + 4) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_D, MID_HC, MID_HWd)
            ijob = ijob + 1
            return
        end if
        ! SWd = D'*SC
        if (ijob == SID_ITER + 5) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_D, MID_SC, MID_SWd)
            ijob = ijob + 1
            return
        end if

        ! - calculate HWdd = D'*H*D = D'*(H*D)
        ! -- HG = H*D
        if (ijob == SID_ITER + 6) then
            call rci_h_multi(iS, task, 'N', m, n, MID_D, MID_HG)
            ijob = ijob + 1
            return
        end if
        ! -- HWdd = D'*HG
        if (ijob == SID_ITER + 7) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_D, MID_HG, MID_HWdd)
            ijob = ijob + 1
            return
        end if

        !  SWdd=D'*S*D
        ! - calculate SWdd = D'*S*D
        ! -- SG = S*D
        if (ijob == SID_ITER + 8) then
            call rci_s_multi(iS, task, 'N', m, n, MID_D, MID_SG)
            ijob = ijob + 1
            return
        end if
        ! -- SWdd = D'*SG
        if (ijob == SID_ITER + 9) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_D, MID_SG, MID_SWdd)
            ijob = SID_COEFF
            return
        end if

        ! using the coeffs. calculated anlytically, we can find
        ! the minimum of the functional in the search direction,
        ! and calculate the energy at that minimum
        if (ijob == SID_UPDATE) then
            call rci_null(iS, task)
            call omm_solve_quartic(coeff(0:4), x_min, ls_fail)
            ! in certain regions of the coeffs. space the line search
            ! gives no minimum--this occurs when there are positive
            ! eigenvalues in the eigenspecturm which are significantly
            ! occupied by our coeffs.  matrix; the only known cure,
            ! unfortunately, is to scale down the entire matrix, thus
            ! returning to a safe region of the coeffs. space.
            if (ls_fail) then
                ijob = SID_LS_FAIL
                e_min_save = 3.0*e_min_save
            else
                ! if the line search is successful, move to the minimum
                e_min = coeff(4)*x_min**4 + &
                        coeff(3)*x_min**3 + &
                        coeff(2)*x_min**2 + &
                        coeff(1)*x_min + &
                        coeff(0)
                e_min_save = e_min
                ls_conv = .true.
                e_diff = 2.0_r8*abs((e_min - e_min_old)/(e_min + e_min_old))
                write (*, '("E_diff: ",E15.3)') e_diff
                if (e_diff <= cg_tol_internal) then
                    conv = .true.
                    call rci_converge(iS, task, conv)
                    return
                end if
                icg = icg + 1
                if (icg > max_iter) then
                    call rci_converge(iS, task, conv)
                    return
                end if
                ijob = SID_UPDATE + 1
            end if
            return
        end if

        ! - C_min = x_min*D + C_min
        if (ijob == SID_UPDATE + 1) then
            call rci_axpy(iS, task, m, n, x_min, MID_D, MID_C)
            ijob = ijob + 1
            return
        end if

        ! - HWdT = HWd'
        if (ijob == SID_UPDATE + 2) then
            call rci_copy(iS, task, 'C', n, n, MID_HWd, MID_HWdT)
            ijob = ijob + 1
            return
        end if

        ! - HW = x_min*HWdT + HW
        if (ijob == SID_UPDATE + 3) then
            call rci_axpy(iS, task, n, n, x_min, MID_HWdT, MID_HW)
            ijob = ijob + 1
            return
        end if

        ! - HW = x_min*HWd + HW
        if (ijob == SID_UPDATE + 4) then
            call rci_axpy(iS, task, n, n, x_min, MID_HWd, MID_HW)
            ijob = ijob + 1
            return
        end if

        ! - HW = x_min^2*HWdd + HW
        if (ijob == SID_UPDATE + 5) then
            call rci_axpy(iS, task, n, n, x_min**2, MID_HWdd, MID_HW)
            ijob = ijob + 1
            return
        end if

        ! - SWdT = SWd'
        if (ijob == SID_UPDATE + 6) then
            call rci_copy(iS, task, 'C', n, n, MID_SWd, MID_SWdT)
            ijob = ijob + 1
            return
        end if

        ! - SW = x_min*SWdT + SW
        if (ijob == SID_UPDATE + 7) then
            call rci_axpy(iS, task, n, n, x_min, MID_SWdT, MID_SW)
            ijob = ijob + 1
            return
        end if

        ! - SW = x_min*SWd + SW
        if (ijob == SID_UPDATE + 8) then
            call rci_axpy(iS, task, n, n, x_min, MID_SWd, MID_SW)
            ijob = ijob + 1
            return
        end if

        ! - SW = x_min^2*SWdd + SW
        if (ijob == SID_UPDATE + 9) then
            call rci_axpy(iS, task, n, n, x_min**2, MID_SWdd, MID_SW)
            ijob = ijob + 1
            return
        end if

        ! recalculate G
        ! - HC = x_min*HG + HC
        if (ijob == SID_UPDATE + 10) then
            call rci_axpy(iS, task, m, n, x_min, MID_HG, MID_HC)
            ijob = ijob + 1
            return
        end if

        ! - SC = x_min*SG + SC
        if (ijob == SID_UPDATE + 11) then
            call rci_axpy(iS, task, m, n, x_min, MID_SG, MID_SC)
            ijob = ijob + 1
            return
        end if

        ! - G = -2 HC*SW
        if (ijob == SID_UPDATE + 12) then
            call rci_gemm(iS, task, 'N', 'N', m, n, n, -2.0_r8, 0.0_r8, &
                          MID_HC, MID_SW, MID_G)
            ijob = ijob + 1
            return
        end if
        ! - G = -2 SC*HW + G
        if (ijob == SID_UPDATE + 13) then
            call rci_gemm(iS, task, 'N', 'N', m, n, n, -2.0_r8, 1.0_r8, &
                          MID_SC, MID_HW, MID_G)
            ijob = ijob + 1
            return
        end if
        ! - G = 4*HC + G
        if (ijob == SID_UPDATE + 14) then
            call rci_axpy(iS, task, m, n, 4.0_r8, MID_HC, MID_G)
            ijob = ijob + 1
            return
        end if

        ! - PG = P*G
        if (ijob == SID_UPDATE + 15) then
            call rci_p_multi(iS, task, 'N', m, n, MID_G, MID_PG)
            ijob = SID_LS_CONV
            return
        end if

        ! calculate lambda
        if (ijob == SID_LS_CONV) then
            call rci_null(iS, task)
            if (ls_conv) then
                ijob = ijob + 1
                return
            else
                lambda = 0.0_r8
                ijob = SID_ITER
                return
            end if
        end if
        ! Work = G
        if (ijob == SID_LS_CONV + 1) then
            call rci_copy(iS, task, 'N', m, n, MID_G, MID_WORK)
            ijob = ijob + 1
            return
        end if
        ! Work = - Gp + Work
        if (ijob == SID_LS_CONV + 2) then
            call rci_axpy(iS, task, m, n, -1.0_r8, MID_Gp, MID_WORK)
            ijob = ijob + 1
            return
        end if
        ! lambda_n = trace(PG*Work)
        if (ijob == SID_LS_CONV + 3) then
            call rci_dot(iS, task, m, n, MID_PG, MID_WORK)
            ijob = ijob + 1
            return
        end if
        ! lambda_d = trace(PGp*Gp)
        if (ijob == SID_LS_CONV + 4) then
            lambda_n = result_in(1)
            call rci_dot(iS, task, m, n, MID_PGp, MID_Gp)
            ijob = ijob + 1
            return
        end if
        ! lambda = lambda_n/lambda_d
        if (ijob == SID_LS_CONV + 5) then
            lambda_d = result_in(1)
            lambda = lambda_n/lambda_d
            call rci_null(iS, task)
            ijob = SID_ITER
            return
        end if

        ! ls_fail
        ! C = 0.5*C
        if (ijob == SID_LS_FAIL) then
            call rci_scale(iS, task, n, n, 0.5_r8, MID_C)
            ijob = ijob + 1
            return
        end if
        ! HW = 0.25*HW
        if (ijob == SID_LS_FAIL + 1) then
            call rci_scale(iS, task, n, n, 0.25_r8, MID_HW)
            ijob = ijob + 1
            return
        end if
        ! SW = 0.25*SW
        if (ijob == SID_LS_FAIL + 2) then
            call rci_scale(iS, task, n, n, 0.25_r8, MID_SW)
            ijob = ijob + 1
            return
        end if
        ! HC = 0.5*HC
        if (ijob == SID_LS_FAIL + 3) then
            call rci_scale(iS, task, m, n, 0.5_r8, MID_HC)
            ijob = ijob + 1
            return
        end if
        ! SC = 0.5*SC
        if (ijob == SID_LS_FAIL + 4) then
            call rci_scale(iS, task, m, n, 0.5_r8, MID_SC)
            ijob = ijob + 1
            return
        end if
        ! G_p = G
        if (ijob == SID_LS_FAIL + 5) then
            call rci_copy(iS, task, 'N', m, n, MID_G, MID_Gp)
            ijob = ijob + 1
            return
        end if
        ! HC = 1.5*G+HC
        if (ijob == SID_LS_FAIL + 6) then
            call rci_axpy(iS, task, m, n, 1.5_r8, MID_G, MID_HC)
            ijob = ijob + 1
            return
        end if
        ! PG = P*G
        if (ijob == SID_LS_FAIL + 7) then
            call rci_p_multi(iS, task, 'N', m, n, MID_G, MID_PG)
            ijob = SID_ITER
            return
        end if

    end subroutine rci_omm

end module ELSI_RCI_OMM
