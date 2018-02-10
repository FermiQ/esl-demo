program esl_demo
 implicit none

 character(len=100) :: input_file

 !Init MPI

 !Init data basis strucutes 

 !Read Input variables and init corresponding data structures
 input_file="sample.inp"
 if (command_argument_count() == 1 ) Then
   call get_command_argument(1, input_file)
 end if

 ! To be able to use YAML
 call f_lib_initialize()
 
 call main(trim(input_file))
 
 !Outputs

 !Release memory
 ! which is not released in final procedure for different types
 call f_lib_finalize()

 !End of the calculation

contains

  subroutine main(input_file)
    use fdf, only : fdf_init, fdf_shutdown, fdf_get
    use esl_hamiltonian_m, only : hamiltonian_t
    use scf_esl, only : scf_t, scf_loop
    use states_esl, only : states_t
    use system_esl, only : system_t
    use numeric_esl, only : init_random
    use smear_esl, only : smear_t
    use elsi_wrapper_esl, only : elsi_t
    use yaml_output
    use esl_next_step_m, only: next_step_setup
    
    character(len = *), intent(in) :: input_file
    type(hamiltonian_t) :: hamiltonian
    type(system_t)      :: system
    type(scf_t)         :: scf
    type(smear_t)       :: smear
    type(states_t)      :: states
    type(elsi_t)        :: elsi
    character(len=100)  :: echo_file,output_file
    integer :: of

    !--------------TEMP --------------------
    integer :: nstates, nspin, nel
    !--------------------------------------

    ! Steps
    integer :: istep, nstep

    echo_file=trim(input_file)//".echo"
    call yaml_map("Input file", trim(input_file))
    call yaml_map("Input file (echo)", trim(echo_file))

    !Parsing the input file
    call fdf_init(input_file, echo_file)

    output_file = fdf_get('output', 'sample.out')
    call yaml_map("Output file", trim(output_file))
    open(newunit=of,file=trim(output_file),action="write")
    call init_random()
    call system%init()
    
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
       
       ! Initialize a new step
       ! This initializes all variables that
       call next_step_setup(system)
       
       call hamiltonian%init(system, states)
       call scf%init()
       call smear%init()

       !SCF loop
       call scf_loop(scf, elsi, hamiltonian, system, states, smear)

    end do

    ! The overhead of having FDF up while running is negligeble.
    ! It serves a good dictionary for parameters.
    call fdf_shutdown()

    close(of)

  end subroutine main
  
end program esl_demo
