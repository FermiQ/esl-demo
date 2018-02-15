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
module ELSI_RCI_DAVIDSON

    use ELSI_RCI_BLAS
    use ELSI_RCI_CONSTANTS
    use ELSI_RCI_DATATYPE
    use ELSI_RCI_PRECISION, only:i4

    implicit none

    public :: rci_davidson

    ! Matrix ID of size m by n
    !&<
    integer(kind=i4), parameter :: MID_Psi     = 1
    integer(kind=i4), parameter :: MID_PsiExt  = 2
    integer(kind=i4), parameter :: MID_HPsi    = 3
    integer(kind=i4), parameter :: MID_HPsiExt = 4
    integer(kind=i4), parameter :: MID_SPsi    = 5
    integer(kind=i4), parameter :: MID_SPsiExt = 6
    !&>

    ! Matrix ID of size n by n
    !&<
    integer(kind=i4), parameter :: MID_HR     = 21
    integer(kind=i4), parameter :: MID_HR12   = 22
    integer(kind=i4), parameter :: MID_HR21   = 23
    integer(kind=i4), parameter :: MID_HR22   = 24
    integer(kind=i4), parameter :: MID_HR1    = 25
    integer(kind=i4), parameter :: MID_HR2    = 26
    integer(kind=i4), parameter :: MID_SR     = 27
    integer(kind=i4), parameter :: MID_SR12   = 28
    integer(kind=i4), parameter :: MID_SR21   = 29
    integer(kind=i4), parameter :: MID_SR22   = 30
    integer(kind=i4), parameter :: MID_SR1    = 31
    integer(kind=i4), parameter :: MID_SR2    = 32
    integer(kind=i4), parameter :: MID_VR     = 33
    integer(kind=i4), parameter :: MID_VRsub  = 34
    !&>

    ! Stage ID
    !&<
    integer(kind=i4), parameter :: SID_INIT   = 0
    integer(kind=i4), parameter :: SID_ITER   = 100
    integer(kind=i4), parameter :: SID_FINISH = 200
    !&>

contains

    subroutine rci_davidson(ijob, iS, task, result_in, m, n_in, e_min, &
                            cvg_tol, max_iter, max_n, long_out)

        implicit none

        !**** INPUT ***********************************!

        logical, intent(in) :: long_out ! print detailed output?

        integer, intent(in) :: m ! size of basis
        integer, intent(in) :: n_in ! number of occupied states
        integer, intent(in) :: max_iter ! max number of iterations
        integer, intent(in) :: max_n ! max number of basis in reduced space

        real(r8), intent(in) :: cvg_tol ! convergence tolerance
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
        integer :: it, iti
        integer, save :: niter ! CG step num.
        integer, save :: next ! number of extended states
        integer, save :: n ! number of current states

        logical, save, allocatable :: Mat_Conv(:)

        real(r8), save, allocatable :: Mat_EW(:)
        real(r8), save, allocatable :: Mat_EWsub(:)
        real(r8), save, allocatable :: Mat_EWold(:)
        real(r8), allocatable :: tmp_vec(:)
        !**********************************************!

        ! -- Init: parameter setup
        if (ijob <= SID_INIT) then
            ijob = SID_INIT + 1

            n = n_in
            niter = 0
            next = 0

            ! n = nband
            call rci_null(iS, task)
            return
        end if

        ! -- HPsi = H*Psi
        if (ijob <= SID_INIT + 1) then
            call rci_h_multi(iS, task, 'N', m, n, MID_Psi, MID_HPsi)
            ijob = ijob + 1
            return
        end if

        ! -- SPsi = S*Psi
        if (ijob == SID_INIT + 2) then
            call rci_s_multi(iS, task, 'N', m, n, MID_Psi, MID_SPsi)
            ijob = ijob + 1
            return
        end if

        ! -- HR = Psi'*HPsi
        if (ijob == SID_INIT + 3) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_Psi, MID_HPsi, MID_HR)
            ijob = ijob + 1
            return
        end if

        ! -- SR = Psi'*SPsi
        if (ijob == SID_INIT + 4) then
            call rci_gemm(iS, task, 'C', 'N', n, n, m, 1.0_r8, 0.0_r8, &
                          MID_Psi, MID_SPsi, MID_SR)
            ijob = ijob + 1
            return
        end if

        ! -- HR * VR = SR * VR * Diag(EW)
        ! - calculate reduced generalized eigenvalue problem
        if (ijob == SID_INIT + 5) then
            call rci_hegv(iS, task, 'V', 'L', n, MID_HR, MID_SR, MID_VR)
            ijob = ijob + 1
            return
        end if

        ! -store old ev for convergence checking
        if (ijob == SID_INIT + 6) then
            Mat_EW = result_in
            next = n
            Mat_EWold = Mat_EW
            call rci_null(iS, task)
            ijob = SID_ITER
            return
        end if

        !-------------------------------------------------
        ! This is the main loop of the Davidson algorithm.

        ! - Redistribution of the vectors
        ! VRsub = VR(Mat_Conv)
        if (ijob == SID_ITER) then
            niter = niter + 1
            call rci_subcol(iS, task, m, next, MID_VR, Mat_Conv, MID_VRsub)
            ijob = ijob + 1
            return
        end if

        ! EWsub = EW(notcvg)
        if (ijob == SID_ITER + 1) then
            deallocate (Mat_EWsub)
            allocate (Mat_EWsub(count(Mat_Conv)))
            iti = 0
            do it = 1, size(Mat_Conv)
                if (.NOT. Mat_Conv(it)) then
                    iti = iti + 1
                    Mat_EWsub(iti) = Mat_EW(it)
                end if
            end do
            call rci_null(iS, task)
            ijob = ijob + 1
            return
        end if

        ! - PsiExt = SPsi * VRsub
        if (ijob == SID_ITER + 2) then
            call rci_gemm(iS, task, 'N', 'N', m, next, n, 1.0_r8, 0.0_r8, &
                          MID_SPsi, MID_VRsub, MID_PsiExt)
            ijob = ijob + 1
            return
        end if

        ! - PsiExt = PsiExt * Diag(EWsub)
        if (ijob == SID_ITER + 3) then
            call rci_colscale(iS, task, m, next, Mat_EWsub, MID_PsiExt)
            ijob = ijob + 1
            return
        end if

        ! - PsiExt = HPsi * VRsub - PsiExt
        if (ijob == SID_ITER + 4) then
            call rci_gemm(iS, task, 'N', 'N', m, next, n, 1.0_r8, -1.0_r8, &
                          MID_HPsi, MID_VRsub, MID_PsiExt)
            ijob = ijob + 1
            return
        end if

        ! - approx solve(PsiExt)
        ! This is a special function only exist for Davidson method
        if (ijob == SID_ITER + 5) then
            call rci_approx_solve(iS, task, m, next, Mat_EWsub, MID_PsiExt)
            ijob = ijob + 1
            return
        end if

        ! - norm(PsiExt)
        if (ijob == SID_ITER + 6) then
            call rci_col_norm(iS, task, m, next, MID_PsiExt)
            ijob = ijob + 1
            return
        end if

        ! - normize(PsiExt)
        if (ijob == SID_ITER + 7) then
            tmp_vec = result_in
            do it = 1, next
                tmp_vec(it) = 1.0_r8/result_in(it)
            end do
            call rci_colscale(iS, task, m, next, tmp_vec, MID_PsiExt)
            ijob = ijob + 1
            return
        end if

        ! -- HPsiExt = H * PsiExt
        if (ijob <= SID_ITER + 8) then
            call rci_h_multi(iS, task, 'N', m, next, MID_PsiExt, &
                             MID_HPsiExt)
            ijob = ijob + 1
            return
        end if

        ! -- HR12 = Psi' * HPsiExt
        if (ijob <= SID_ITER + 9) then
            call rci_gemm(iS, task, 'T', 'N', n, next, m, 1.0_r8, 0.0_r8, &
                          MID_Psi, MID_HPsiExt, MID_HR12)
            ijob = ijob + 1
            return
        end if

        ! -- HR22 = PsiExt' * HPsiExt
        if (ijob <= SID_ITER + 10) then
            call rci_gemm(iS, task, 'T', 'N', next, next, m, &
                          1.0_r8, 0.0_r8, &
                          MID_PsiExt, MID_HPsiExt, MID_HR22)
            ijob = ijob + 1
            return
        end if

        ! -- HR1 = [HR HR12]
        if (ijob == SID_ITER + 11) then
            call rci_cbind(iS, task, n, n, next, MID_HR, MID_HR12, MID_HR1)
            ijob = ijob + 1
            return
        end if

        ! -- HR2 = [HR21 HR22]
        if (ijob == SID_ITER + 12) then
            call rci_cbind(iS, task, next, n, next, &
                           MID_HR21, MID_HR22, MID_HR2)
            ijob = ijob + 1
            return
        end if

        ! -- HR = [HR1; HR2]
        if (ijob == SID_ITER + 13) then
            call rci_rbind(iS, task, n, next, n + next, &
                           MID_HR1, MID_HR2, MID_HR)
            ijob = ijob + 1
            return
        end if

        ! -- SPsiExt = S*PsiExt
        if (ijob == SID_ITER + 14) then
            call rci_s_multi(iS, task, 'N', m, next, &
                             MID_PsiExt, MID_SPsiExt)
            ijob = ijob + 1
            return
        end if

        ! -- SR12 = Psi' * SPsiExt
        if (ijob <= SID_ITER + 15) then
            call rci_gemm(iS, task, 'T', 'N', n, next, m, 1.0_r8, 0.0_r8, &
                          MID_Psi, MID_SPsiExt, MID_SR12)
            ijob = ijob + 1
            return
        end if

        ! -- SR22 = PsiExt' * SPsiExt
        if (ijob <= SID_ITER + 16) then
            call rci_gemm(iS, task, 'T', 'N', next, next, m, &
                          1.0_r8, 0.0_r8, &
                          MID_PsiExt, MID_SPsiExt, MID_SR22)
            ijob = ijob + 1
            return
        end if

        ! -- SR21 = SR12^(T)
        if (ijob == SID_ITER + 17) then
            call rci_copy(iS, task, 'C', next, n, MID_SR12, MID_SR21)
            ijob = ijob + 1
            return
        end if

        ! -- SR1 = [SR SR12]
        if (ijob == SID_ITER + 18) then
            call rci_cbind(iS, task, n, n, next, MID_SR, MID_SR12, MID_SR1)
            ijob = ijob + 1
            return
        end if

        ! -- SR2 = [SR21 SR22]
        if (ijob == SID_ITER + 19) then
            call rci_cbind(iS, task, next, n, next, &
                           MID_SR21, MID_SR22, MID_SR2)
            ijob = ijob + 1
            return
        end if

        ! -- SR = [SR1; SR2]
        if (ijob == SID_ITER + 20) then
            call rci_rbind(iS, task, n, next, n + next, &
                           MID_SR1, MID_SR2, MID_SR)
            ijob = ijob + 1
            return
        end if

        ! -- HR * VR = SR * VR * Diag(EW)
        ! - calculate reduced generalized eigenvalue problem
        if (ijob == SID_ITER + 21) then
            call rci_hegv(iS, task, 'V', 'L', n + next, &
                          MID_HR, MID_SR, MID_VR)
            ijob = ijob + 1
            return
        end if

        ! -store old ev for convergence checking
        if (ijob == SID_ITER + 22) then
            Mat_EW = result_in
            Mat_Conv = (abs(Mat_EW(1:n) - Mat_EWold(1:n)) < cvg_tol)
            next = count(.NOT. Mat_Conv)
            n = n + next
            Mat_EWold = Mat_EW
            call rci_null(iS, task)
            if (next == 0 .OR. n > max_n .OR. niter == max_iter) then
                ijob = SID_FINISH
                return
            else
                ijob = SID_ITER
                return
            end if
        end if

        ! - Finish the iteration
        ! - PsiExt = Psi * VR
        if (ijob == SID_FINISH) then
            call rci_gemm(iS, task, 'N', 'N', m, n, n, 1.0_r8, 0.0_r8, &
                          MID_Psi, MID_VR, MID_PsiExt)
            ijob = ijob + 1
            return
        end if

        ! - Psi = PsiExt
        if (ijob == SID_FINISH + 1) then
            call rci_copy(iS, task, 'N', m, n, MID_PsiExt, MID_Psi)
            ijob = ijob + 1
            return
        end if

        ! - Converge or Stop
        if (ijob == SID_FINISH + 1) then
            call rci_converge(iS, task, next == 0)
            return
        end if

    end subroutine rci_davidson

end module ELSI_RCI_DAVIDSON
