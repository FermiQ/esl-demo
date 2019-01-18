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
!! This module contains constants used in ELSI RCI.
!!
module ELSI_RCI_CONSTANTS

    use ELSI_RCI_PRECISION, only: i4

    implicit none

    ! Algorithms
    integer(kind=i4), parameter :: ELSI_RCI_JD = 0
    integer(kind=i4), parameter :: ELSI_RCI_CG = 1
    ! ...

    integer(kind=i4), parameter :: ELSI_RCI_NULL     = 0
    integer(kind=i4), parameter :: ELSI_RCI_CONVERGE = 1
    integer(kind=i4), parameter :: ELSI_RCI_STOP     = 12
    integer(kind=i4), parameter :: ELSI_RCI_H_MULTI  = 2
    integer(kind=i4), parameter :: ELSI_RCI_S_MULTI  = 3
    integer(kind=i4), parameter :: ELSI_RCI_P_MULTI  = 4
    integer(kind=i4), parameter :: ELSI_RCI_GEMM     = 5
    integer(kind=i4), parameter :: ELSI_RCI_AXPY     = 6
    integer(kind=i4), parameter :: ELSI_RCI_COPY     = 7
    integer(kind=i4), parameter :: ELSI_RCI_TRACE    = 8
    integer(kind=i4), parameter :: ELSI_RCI_DOT      = 9
    integer(kind=i4), parameter :: ELSI_RCI_SCALE    = 10
    integer(kind=i4), parameter :: ELSI_RCI_COLSCALE = 15
    integer(kind=i4), parameter :: ELSI_RCI_ROWSCALE = 16
    integer(kind=i4), parameter :: ELSI_RCI_HEGV     = 11
    integer(kind=i4), parameter :: ELSI_RCI_CBIND    = 13
    integer(kind=i4), parameter :: ELSI_RCI_RBIND    = 14
    integer(kind=i4), parameter :: ELSI_RCI_COL_NORM = 17
    integer(kind=i4), parameter :: ELSI_RCI_SUBCOL   = 18

    integer(kind=i4), parameter :: ELSI_RCI_APPROX_SOLVE = 31

end module ELSI_RCI_CONSTANTS
