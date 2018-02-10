program esl_demo
  use yaml_output

#ifdef WITH_LUA
  use flook, only : luaState
  use esl_flook_if_m
#endif

  implicit none

  ! Input-file
  character(len=256) :: input_file
  integer :: out_unit
  
#ifdef WITH_LUA
  ! LUA-handle, we really do need it up here to be able to use it
  type(luaState) :: LUA
#endif


  ! Initialization of the program
  call esl_init()

  ! Call ESL
  call esl_main()

  ! Final end of stuff, fdf, flib, MPI, etc.
  call esl_end()

contains

  subroutine esl_main()
    use esl_hamiltonian_m, only : hamiltonian_t
    use esl_scf_m, only : scf_t, scf_loop
    use esl_states_m, only : states_t
    use esl_system_m, only : system_t
    use esl_smear_m, only : smear_t
    use elsi_wrapper_esl, only : elsi_t
    use esl_next_step_m, only: next_step_setup

    type(hamiltonian_t) :: hamiltonian
    type(system_t)      :: system
    type(scf_t)         :: scf
    type(smear_t)       :: smear
    type(states_t)      :: states
    type(elsi_t)        :: elsi

    !--------------TEMP --------------------
    integer :: nstates, nspin, nel
    !--------------------------------------

    ! Steps
    integer :: istep, nstep

    call system%init()

#ifdef WITH_LUA
    ! Call lua at given state
    call flook_if_call(LUA, LUA_INITIALIZE)
#endif

    !--------------TEMP --------------------
    nstates = 1
    nspin = 1
    nel = 1
    !---------------------------------------
    call states%init(system%basis, nstates, nspin, 1, nel)
    call states%summary()

    ! TODO Nstep should probably be read from another entity
    ! Or even better be a logical (do while ( step_finished )
    nstep = 1
    do istep = 1, nstep ! steps to run through for the SCF
      
      ! this step is not adhearing to anything, it could be
      ! md-steps, or whatever.

#ifdef WITH_LUA
      ! Call lua just before we are to initialize a new step.
      call flook_if_call(LUA, LUA_INIT_STEP)
#endif

      ! This initializes all variables that
      call next_step_setup(system)

      call hamiltonian%init(system, states)
      call scf%init()
      call smear%init()

      !SCF loop
      call scf_loop(scf, elsi, hamiltonian, system, states, smear)

#ifdef WITH_LUA
      ! Call lua after the forces has been calculated.
      call flook_if_call(LUA, LUA_FORCES)
#endif

      ! Insert molecular dynamics, etc.
      
#ifdef WITH_LUA
      ! Call lua after the SCF-loop.
      ! At this point we can change variables, i.e. perform MD etc.
      call flook_if_call(LUA, LUA_NEXT_STEP)
#endif

    end do

  end subroutine esl_main

  subroutine esl_init()
    use fdf, only : fdf_init, fdf_get
    use esl_numeric_m, only : init_random

    character(len=256) :: echo_file

    !Init MPI

    !Init data basis strucutes 

    !Read Input variables and init corresponding data structures
    input_file = "sample.inp"
    if ( command_argument_count() == 1 ) Then
      call get_command_argument(1, input_file)
    end if

    ! To be able to use YAML
    call f_lib_initialize()

    ! TODO abstraction with MPI (all nodes will write to the same file)
    echo_file = trim(input_file)//".echo"
    call yaml_map("Input file", trim(input_file))
    call yaml_map("Input file (echo)", trim(echo_file))

    ! Open fdf and parse input file
    call fdf_init(input_file, echo_file)

    ! Read the output file
    output_file = fdf_get('Output', 'sample.out')
    call yaml_map("Output file", trim(output_file))

    open(newunit=out_unit, file=trim(output_file), action="write")

    ! Initialize random sequences
    call init_random()

#ifdef WITH_LUA
    ! Initialize Lua interface
    call flook_if_init(LUA)
#endif

  end subroutine esl_init

  subroutine esl_end()
    use fdf, only : fdf_shutdown

    !Release memory
    ! which is not released in final procedure for different types
    call f_lib_finalize()

    ! The overhead of having FDF up while running is negligeble.
    ! It serves a good dictionary for parameters.
    call fdf_shutdown()

    close(out_unit)

  end subroutine esl_end

end program esl_demo
