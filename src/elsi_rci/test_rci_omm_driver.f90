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
!! This program tests elsi_dm_real.
!!
program test_dm_real

    use ELSI
    use ELSI_RCI
    use ELSI_RCI_CONSTANTS
    use ELSI_RCI_OMM
    use ELSI_RCI_PRECISION, only:r8, i4
    use MatrixSwitch

    implicit none

    include "mpif.h"

    character(128) :: arg1 ! solver
    character(128) :: arg2 ! H file
    character(128) :: arg3 ! S file
    character(128) :: arg4 ! make check?

    character(5)   :: m_storage
    character(3)   :: m_operation

    integer(kind=i4) :: n_proc
    integer(kind=i4) :: nprow
    integer(kind=i4) :: npcol
    integer(kind=i4) :: myid
    integer(kind=i4) :: mpi_comm_global
    integer(kind=i4) :: mpierr
    integer(kind=i4) :: blk
    integer(kind=i4) :: blacs_ctxt
    integer(kind=i4) :: n_states
    integer(kind=i4) :: matrix_size
    integer(kind=i4) :: solver
    integer(kind=i4) :: task
    integer(kind=i4) :: i
    integer(kind=i4) :: j
    integer(kind=i4) :: m
    integer(kind=i4) :: n
    integer(kind=i4) :: l_rows
    integer(kind=i4) :: l_cols
    integer(kind=i4) :: ijob
    integer(kind=i4) :: max_iter = 10000

    real(kind=r8) :: n_electrons
    real(kind=r8) :: e_test = 0.0_r8
    real(kind=r8) :: e_ref = 0.0_r8
    real(kind=r8) :: e_tol = 0.0_r8
    real(kind=r8) :: t1
    real(kind=r8) :: t2
    real(kind=r8) :: rnum
    real(kind=r8) :: e_min = 0.0_r8
    real(kind=r8) :: cg_tol = 0.0001_r8

    logical :: make_check = .false. ! Are we running "make check"?
    logical :: long_out = .false.

    real(kind=r8), allocatable :: ham(:, :)
    real(kind=r8), allocatable :: ovlp(:, :)
    real(kind=r8), allocatable :: result_in(:)

    type(elsi_handle)    :: e_h
    type(elsi_rw_handle) :: rw_h
    type(rci_instr)      :: iS

    type(matrix) :: H
    type(matrix) :: S
    type(matrix), allocatable :: Work(:)

    ! VY: Reference values from calculations on August 31, 2017.
    real(kind=r8), parameter :: e_omm = -1833.07932666692_r8

    ! Initialize MPI
    call MPI_Init(mpierr)
    mpi_comm_global = MPI_COMM_WORLD
    call MPI_Comm_size(mpi_comm_global, n_proc, mpierr)
    call MPI_Comm_rank(mpi_comm_global, myid, mpierr)

    ! Read command line arguments
    if (COMMAND_ARGUMENT_COUNT() == 3) then
        call GET_COMMAND_ARGUMENT(1, arg1)
        call GET_COMMAND_ARGUMENT(2, arg2)
        call GET_COMMAND_ARGUMENT(3, arg3)
        read (arg1, *) solver
    else
        if (myid == 0) then
            write (*, '("  #############################################")')
            write (*, '("  ## Wrong number of command line arguments!!##")')
            write (*, '("  ## Arg#1: Choice of solver.                ##")')
            write (*, '("  ##        (OMM_RCI = 1)                    ##")')
            write (*, '("  ## Arg#2: H matrix file.                   ##")')
            write (*, '("  ## Arg#3: S matrix file.                   ##")')
            write (*, '("  #############################################")')
            call MPI_Abort(mpi_comm_global, 0, mpierr)
            stop
        endif
    endif

    if (myid == 0) then
        e_tol = 1.0e-8_r8
        write (*, '("  ####################################")')
        write (*, '("  ##     ELSI RCI TEST PROGRAMS     ##")')
        write (*, '("  ####################################")')
        write (*, *)
        if (solver == 2) then
            write (*, '("  This test program performs the following computational steps:")')
            write (*, *)
            write (*, '("  1) Reads Hamiltonian and overlap matrices;")')
            write (*, '("  2) Computes the density matrix with OMM_RCI interface.")')
            write (*, *)
            write (*, '("  Now start testing  elsi_dm_real + OMM_RCI")')
            e_ref = e_omm
        endif
        write (*, *)
    endif

    ! Set up square-like processor grid
    do npcol = nint(sqrt(real(n_proc))), 2, -1
        if (mod(n_proc, npcol) == 0) exit
    enddo
    nprow = n_proc/npcol

    ! Set block size
    blk = 32

    ! Set up BLACS
    blacs_ctxt = mpi_comm_global
    call BLACS_Gridinit(blacs_ctxt, 'r', nprow, npcol)

    ! Read H and S matrices
    call elsi_init_rw(rw_h, 0, 1, 0, 0.0_r8)
    call elsi_set_rw_mpi(rw_h, mpi_comm_global)
    call elsi_set_rw_blacs(rw_h, blacs_ctxt, blk)

    call elsi_read_mat_dim(rw_h, arg2, n_electrons, matrix_size, &
                           l_rows, l_cols)

    allocate (ham(l_rows, l_cols))
    allocate (ovlp(l_rows, l_cols))

    t1 = MPI_Wtime()

    call elsi_read_mat_real(rw_h, arg2, ham)
    call elsi_read_mat_real(rw_h, arg3, ovlp)

    call elsi_finalize_rw(rw_h)

    t2 = MPI_Wtime()

    if (myid == 0) then
        write (*, '("  Finished reading H and S matrices")')
        write (*, '("  | Time :",F10.3,"s")') t2 - t1
        write (*, *)
    endif

    ! Initialize ELSI
    n_states = int(n_electrons, kind=i4)

    call elsi_init(e_h, solver, 1, 0, matrix_size, n_electrons, n_states)
    call elsi_set_mpi(e_h, mpi_comm_global)
    call elsi_set_blacs(e_h, blacs_ctxt, blk)

    call m_register_pdbc(H, ham, e_h%sc_desc)
    call m_register_pdbc(S, ovlp, e_h%sc_desc)

    t1 = MPI_Wtime()

    allocate (Work(28))

    m = matrix_size
    n = n_states

    m_storage = 'pddbc'
    m_operation = 'lap'
    do i = 1, 11
        call m_allocate(Work(i), m, n, m_storage)
    end do
    do i = 21, 28
        call m_allocate(Work(i), n, n, m_storage)
    end do

    ! Initialize wave functions
    do i = 1, m
        do j = 1, n
            call random_number(rnum)
            call m_set_element(Work(1), i, j, rnum, 0.0_r8, m_operation)
        end do
    end do
    call m_scale(Work(1), 1.0d-2/sqrt(real(m, r8)), m_operation)

    allocate (result_in(1))
    ijob = -1
    do
        call rci_omm(ijob, iS, task, result_in, m, n, &
                     e_min, cg_tol, max_iter, long_out)

        select case (task)
        case (ELSI_RCI_NULL)
        case (ELSI_RCI_CONVERGE)
            exit
        case (ELSI_RCI_H_MULTI)
            call mm_multiply(H, iS%TrH, Work(iS%Aidx), &
                             'n', Work(iS%Bidx), 1.0_r8, 0.0_r8, &
                             m_operation)
        case (ELSI_RCI_S_MULTI)
            call mm_multiply(S, iS%TrS, Work(iS%Aidx), &
                             'n', Work(iS%Bidx), 1.0_r8, 0.0_r8, &
                             m_operation)
        case (ELSI_RCI_P_MULTI)
            ! No preconditioner
            call m_add(Work(iS%Aidx), 'n', Work(iS%Bidx), 1.0_r8, &
                       0.0_r8, m_operation)
        case (ELSI_RCI_GEMM)
            call mm_multiply(Work(iS%Aidx), iS%TrA, Work(iS%Bidx), &
                             iS%TrB, Work(iS%Cidx), iS%alpha, iS%beta, &
                             m_operation)
        case (ELSI_RCI_AXPY)
            call m_add(Work(iS%Aidx), 'n', Work(iS%Bidx), iS%alpha, &
                       1.0_r8, m_operation)
        case (ELSI_RCI_COPY)
            call m_add(Work(iS%Aidx), iS%TrA, Work(iS%Bidx), 1.0_r8, &
                       0.0_r8, m_operation)
        case (ELSI_RCI_TRACE)
            call m_trace(Work(iS%Aidx), result_in(1), m_operation)
        case (ELSI_RCI_DOT)
            call mm_trace(Work(iS%Aidx), Work(iS%Bidx), result_in(1), &
                          m_operation)
        case (ELSI_RCI_SCALE)
            call m_scale(Work(iS%Aidx), iS%alpha, m_operation)
        case default
        end select

    end do

    deallocate (result_in)

    do i = 1, 11
        call m_deallocate(Work(i))
    end do
    do i = 21, 28
        call m_deallocate(Work(i))
    end do

    t2 = MPI_Wtime()

    if (myid == 0) then
        write (*, '("  Finished SCF")')
        write (*, '("  | Time :",F10.3,"s")') t2 - t1
        write (*, *)
        write (*, '("  Finished test program")')
        write (*, *)
        if (make_check) then
            if (abs(e_test - e_ref) < e_tol) then
                write (*, '("  Passed.")')
            else
                write (*, '("  Failed!!")')
            endif
            write (*, *)
        endif
    endif

    ! Finalize ELSI
    call elsi_finalize(e_h)

    deallocate (ham)
    deallocate (ovlp)

    call MPI_Finalize(mpierr)

end program
