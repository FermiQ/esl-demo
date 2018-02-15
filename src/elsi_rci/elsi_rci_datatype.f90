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
!! This module contains variables accessible in ELSI RCI and related modules.
!!
module ELSI_RCI_DATATYPE

    use, intrinsic :: ISO_C_BINDING
    use ELSI_RCI_PRECISION

    implicit none

    private

    !**** TYPES *************************************!
    type, public :: elsi_rci_handle

        integer(kind=i4) :: method
        integer(kind=i4) :: tolerance

        ! ...

    end type

    ! RCI Instruction type
    ! Instructions for H_Multi, S_Multi, P_Multi, GEMM, AXPY, COPY,
    ! TRACE, INNER
    ! H: Hamiltonian; S: Overlapping matrix; P: Preconditioner
    ! A, B, C: General matrices
    type, public :: rci_instr 
        character :: jobz, uplo ! job char; and upper or lower 
        character :: TrH, TrS, TrP, TrA, TrB ! Operation for H, S, P, A, B 

        integer   :: m, n ! size of the output matrix
        integer   :: n1, n2 ! size of matrices in binding
        integer   :: k    ! size of the intermedia multiplication 

        real(kind=r8)  :: alpha, beta ! coefficients
        real(kind=r8),allocatable  :: alphavec(:) ! vector coefficients

        logical,allocatable  :: Idxs(:) ! indicator vector for subcol/row

        integer   :: Aidx, Bidx, Cidx ! indices for matrix A, B and C
    end type

end module ELSI_RCI_DATATYPE
