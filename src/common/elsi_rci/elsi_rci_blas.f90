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

module ELSI_RCI_BLAS

   use ELSI_RCI_CONSTANTS
   use ELSI_RCI_DATATYPE
   use ELSI_RCI_PRECISION, only: r8,i4

   implicit none


contains

    ! Null
    ! - rci_null(iS,task)

    ! Converge flag
    ! - rci_converge(iS,task)

    ! B = H^(trH) * A
    ! - rci_h_multi(iS,task,trH,m,n,Aidx,Bidx)

    ! B = S^(trS) * A
    ! - rci_s_multi(iS,task,trS,m,n,Aidx,Bidx)

    ! B = P^(trP) * A
    ! - rci_p_multi(iS,task,trP,m,n,Aidx,Bidx)

    ! C = alpha * A^(trA) * B^(trB) + beta * C
    ! - rci_gemm(iS,task,trA,trB,m,n,k,alpha,beta,Aidx,Bidx,Cidx)

    ! B = alpha * A + B
    ! - rci_axpy(iS,task,m,n,alpha,Aidx,Bidx)

    ! B = A^(trA)
    ! - rci_copy(iS,task,trA,m,n,Aidx,Bidx)

    ! A = alpha * A
    ! - rci_scale(iS,task,m,n,alpha,Aidx)

    ! res = trace(A)
    ! - rci_trace(iS,task,m,n,Aidx)

    ! res = trace(A * B)
    ! - rci_dot(iS,task,m,n,Aidx,Bidx)

    ! A * C = B * C * Diag(resvec)
    ! - rci_hegv(iS,task,jobz,uplo,n,Aidx,Bidx,Cidx)

    ! C = [A B]
    ! - rci_cbind(iS,task,m,n1,n2,Aidx,Bidx,Cidx)

    ! C = [A; B]
    ! - rci_rbind(iS,task,n1,n2,n,Aidx,Bidx,Cidx)

    ! A = A * Diag(alphavec)
    ! - rci_colscale(iS,task,m,n,alphavec,Aidx)

    ! A = Diag(alphavec) * A
    ! - rci_rowscale(iS,task,m,n,alphavec,Aidx)

    ! resvec = diag( A' * A ) 
    ! - rci_col_norm(iS,task,m,n,Aidx)

    ! B = A(:,Idxs)
    ! - rci_subcol(iS,task,m,n,Aidx,Idxs,Bidx)

    ! Approximately solve (H - lambda * S)^(-1)
    ! - rci_approx_solve(iS,task,m,n,alphavec,Aidx)

subroutine rci_null(iS,task)

    !**** FUNCTION ********************************!
    ! Not do anything

    implicit none

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_NULL

end subroutine rci_null

subroutine rci_converge(iS,task,conv_flag)

    !**** FUNCTION ********************************!
    ! Not do anything

    implicit none

    !**** INPUT ***********************************!
    logical, intent(in) :: conv_flag ! convergence flag

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    if ( conv_flag ) then
        task = ELSI_RCI_CONVERGE
    else
        task = ELSI_RCI_STOP
    end if

end subroutine rci_converge

subroutine rci_h_multi(iS,task,trH,m,n,Aidx,Bidx)

    !**** FUNCTION ********************************!
    ! B = op(H) * A
    ! op() denotes the operation of normal(N), transpose(T), conjugate(C)

    implicit none

    !**** INPUT ***********************************!
    character, intent(in) :: trH ! transpose of H 
    integer, intent(in) :: m ! height of B
    integer, intent(in) :: n ! width of B
    integer, intent(in) :: Aidx, Bidx ! indices for matrix A and B

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_H_MULTI
    iS%TrH   = trH
    iS%m     = m
    iS%n     = n
    iS%Aidx  = Aidx
    iS%Bidx  = Bidx

end subroutine rci_h_multi

subroutine rci_s_multi(iS,task,trS,m,n,Aidx,Bidx)

    !**** FUNCTION ********************************!
    ! B = op(S) * A
    ! op() denotes the operation of normal(N), transpose(T), conjugate(C)

    implicit none

    !**** INPUT ***********************************!
    character, intent(in) :: trS ! transpose of S
    integer, intent(in) :: m ! height of B
    integer, intent(in) :: n ! width of B
    integer, intent(in) :: Aidx, Bidx ! indices for matrix A and B

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_S_MULTI
    iS%TrS   = trS
    iS%m     = m
    iS%n     = n
    iS%Aidx  = Aidx
    iS%Bidx  = Bidx

end subroutine rci_s_multi

subroutine rci_p_multi(iS,task,trP,m,n,Aidx,Bidx)

    !**** FUNCTION ********************************!
    ! B = op(P) * A
    ! op() denotes the operation of normal(N), transpose(T), conjugate(C)

    implicit none

    !**** INPUT ***********************************!
    character, intent(in) :: trP ! transpose of P
    integer, intent(in) :: m ! height of B
    integer, intent(in) :: n ! width of B
    integer, intent(in) :: Aidx, Bidx ! indices for matrix A and B

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_P_MULTI
    iS%TrP   = trP
    iS%m     = m
    iS%n     = n
    iS%Aidx  = Aidx
    iS%Bidx  = Bidx

end subroutine rci_p_multi

subroutine rci_gemm(iS,task,trA,trB,m,n,k,alpha,beta,Aidx,Bidx,Cidx)

    !**** FUNCTION ********************************!
    ! C = alpha * op(A)*op(B) + beta * C
    ! op() denotes the operation of normal(N), transpose(T), conjugate(C)

    implicit none

    !**** INPUT ***********************************!
    character, intent(in) :: trA, trB ! transpose of A and B
    integer, intent(in) :: m ! height of C
    integer, intent(in) :: n ! width of C
    integer, intent(in) :: k ! inner size of op(A)*op(B)
    real(r8), intent(in) :: alpha, beta
    integer, intent(in) :: Aidx, Bidx, Cidx ! indices for matrix A, B and C

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_GEMM
    iS%TrA   = trA
    iS%TrB   = trB
    iS%m     = m
    iS%n     = n
    iS%k     = k
    iS%alpha = alpha
    iS%beta  = beta
    iS%Aidx  = Aidx
    iS%Bidx  = Bidx
    iS%Cidx  = Cidx

end subroutine rci_gemm

subroutine rci_axpy(iS,task,m,n,alpha,Aidx,Bidx)

    !**** FUNCTION ********************************!
    ! B = alpha * A + B

    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m ! height of B
    integer, intent(in) :: n ! width of B
    real(r8), intent(in) :: alpha
    integer, intent(in) :: Aidx, Bidx ! indices for matrix A and B

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_AXPY
    iS%m     = m
    iS%n     = n
    iS%alpha = alpha
    iS%Aidx  = Aidx
    iS%Bidx  = Bidx

end subroutine rci_axpy

subroutine rci_copy(iS,task,trA,m,n,Aidx,Bidx)

    !**** FUNCTION ********************************!
    ! B = op(A)

    implicit none

    !**** INPUT ***********************************!
    character, intent(in) :: trA ! transpose of A
    integer, intent(in) :: m ! height of B
    integer, intent(in) :: n ! width of B
    integer, intent(in) :: Aidx, Bidx ! indices for matrix A and B

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_COPY
    iS%TrA   = trA
    iS%m     = m
    iS%n     = n
    iS%Aidx  = Aidx
    iS%Bidx  = Bidx

end subroutine rci_copy

subroutine rci_scale(iS,task,m,n,alpha,Aidx)

    !**** FUNCTION ********************************!
    ! A = alpha*A

    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m ! height of B
    integer, intent(in) :: n ! width of B
    real(r8), intent(in) :: alpha
    integer, intent(in) :: Aidx ! indices for matrix A

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_SCALE
    iS%m     = m
    iS%n     = n
    iS%alpha = alpha
    iS%Aidx  = Aidx

end subroutine rci_scale

subroutine rci_trace(iS,task,m,n,Aidx)

    !**** FUNCTION ********************************!
    ! res = trace(A)

    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m ! height of A
    integer, intent(in) :: n ! width of A
    integer, intent(in) :: Aidx ! indices for matrix A

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_TRACE
    iS%m     = m
    iS%n     = n
    iS%Aidx  = Aidx

end subroutine rci_trace

subroutine rci_dot(iS,task,m,n,Aidx,Bidx)

    !**** FUNCTION ********************************!
    ! res = trace(A*B) = dot(A(:),B(:))

    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m ! height of A
    integer, intent(in) :: n ! width of A
    integer, intent(in) :: Aidx, Bidx ! indices for matrix A and B

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_DOT
    iS%m     = m
    iS%n     = n
    iS%Aidx  = Aidx
    iS%Bidx  = Bidx

end subroutine rci_dot

subroutine rci_hegv(iS,task,jobz,uplo,n,Aidx,Bidx,Cidx)

    !**** FUNCTION ********************************!
    ! A * C = B * C * Diag(res)

    implicit none

    !**** INPUT ***********************************!
    character, intent(in) :: jobz ! with/without eigenvector
    character, intent(in) :: uplo ! upper or lower part of A
    integer, intent(in) :: n ! size of A
    integer, intent(in) :: Aidx, Bidx, Cidx ! indices for matrix A, B, C

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_HEGV
    iS%jobz  = jobz
    iS%uplo  = uplo
    iS%n     = n
    iS%Aidx  = Aidx
    iS%Bidx  = Bidx
    iS%Cidx  = Cidx

end subroutine rci_hegv

subroutine rci_cbind(iS,task,m,n1,n2,Aidx,Bidx,Cidx)

    !**** FUNCTION ********************************!
    ! C = [A B]

    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m ! height of A, B, C
    integer, intent(in) :: n1, n2 ! width of A and B
    integer, intent(in) :: Aidx, Bidx, Cidx ! indices for matrix A, B, C

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_CBIND
    iS%m  = m
    iS%n  = n1+n2
    iS%n1  = n1
    iS%n2  = n2
    iS%Aidx  = Aidx
    iS%Bidx  = Bidx
    iS%Cidx  = Cidx

end subroutine rci_cbind

subroutine rci_rbind(iS,task,n1,n2,n,Aidx,Bidx,Cidx)

    !**** FUNCTION ********************************!
    ! C = [A; B]

    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: n ! width of A, B, C
    integer, intent(in) :: n1, n2 ! height of A and B
    integer, intent(in) :: Aidx, Bidx, Cidx ! indices for matrix A, B, C

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_RBIND
    iS%m  = n1+n2
    iS%n1  = n1
    iS%n2  = n2
    iS%n  = n
    iS%Aidx  = Aidx
    iS%Bidx  = Bidx
    iS%Cidx  = Cidx

end subroutine rci_rbind

subroutine rci_colscale(iS,task,m,n,alphavec,Aidx)

    !**** FUNCTION ********************************!
    ! A = A * Diag(alpha)

    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m ! height of A
    integer, intent(in) :: n ! width of A
    real(r8), intent(in) :: alphavec(:)
    integer, intent(in) :: Aidx ! indices for matrix A

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_COLSCALE
    iS%m        = m
    iS%n        = n
    iS%alphavec = alphavec
    iS%Aidx     = Aidx

end subroutine rci_colscale

subroutine rci_rowscale(iS,task,m,n,alphavec,Aidx)

    !**** FUNCTION ********************************!
    ! A = Diag(alpha) * A

    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m ! height of A
    integer, intent(in) :: n ! width of A
    real(r8), intent(in) :: alphavec(:)
    integer, intent(in) :: Aidx ! indices for matrix A

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_ROWSCALE
    iS%m        = m
    iS%n        = n
    iS%alphavec = alphavec
    iS%Aidx     = Aidx

end subroutine rci_rowscale

subroutine rci_col_norm(iS,task,m,n,Aidx)

    !**** FUNCTION ********************************!
    ! resvec = diag( A' * A ) 

    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m ! height of A
    integer, intent(in) :: n ! width of A
    integer, intent(in) :: Aidx ! indices for matrix A

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_COL_NORM
    iS%m        = m
    iS%n        = n
    iS%Aidx     = Aidx

end subroutine rci_col_norm

subroutine rci_subcol(iS,task,m,n,Aidx,Idxs,Bidx)

    !**** FUNCTION ********************************!
    ! B = A(:,Idxs)

    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m ! height of B
    integer, intent(in) :: n ! width of B
    integer, intent(in) :: Aidx, Bidx ! indices for matrix A and B
    logical, intent(in) :: Idxs(:) ! indices

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_SUBCOL
    iS%m        = m
    iS%n        = n
    iS%Aidx     = Aidx
    iS%Idxs     = Idxs
    iS%Bidx     = Bidx

end subroutine rci_subcol


!--------------------------------------------------
! Special routines
subroutine rci_approx_solve(iS,task,m,n,alphavec,Aidx)

    !**** FUNCTION ********************************!
    ! A = g_psi( A, alphavec )
    ! This is a special function only exist in Davidson method

    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m ! height of A
    integer, intent(in) :: n ! width of A
    real(r8), intent(in) :: alphavec(:)
    integer, intent(in) :: Aidx ! indices for matrix A

    !**** OUTPUT ***********************************!
    type(rci_instr), intent(out)   :: iS
    integer(kind=i4), intent(out) :: task

    task = ELSI_RCI_APPROX_SOLVE
    iS%m        = m
    iS%n        = n
    iS%alphavec = alphavec
    iS%Aidx     = Aidx

end subroutine rci_approx_solve

end module ELSI_RCI_BLAS
