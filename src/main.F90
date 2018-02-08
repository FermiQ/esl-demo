program esl_demo
 use fdf, only : fdf_init, fdf_shutdown, fdf_string
 use hamiltonian_esl, only : hamiltonian_t
 use iso_fortran_env, only : ou=>OUTPUT_UNIT
 use scf_esl, only : scf_t, scf_loop
 use states_esl
 use system_esl, only : system_t
 use numeric_esl, only : init_random
 use smear_esl, only : smear_t
 use elsi_wrapper_esl, only : elsi_t

 implicit none

 type(hamiltonian_t) :: hamiltonian
 type(system_t)      :: system
 type(scf_t)         :: scf
 type(smear_t)         :: smear
 type(states_t)         :: states
 type(elsi_t)         :: elsi
 character(len=100) :: input_file,echo_file,output_file
 integer :: of

 !--------------TEMP --------------------
 integer :: nstates, nspin
 !--------------------------------------


 !Init MPI

 !Init data basis strucutes 

 !Read Input variables and init corresponding data structures
 input_file="sample.inp"
 if (command_argument_count() == 1 ) Then
   call get_command_argument(1, input_file)
 end if

 echo_file=trim(input_file)//".echo"
 write(ou,'(a)')"reading instructions from: "//trim(input_file)//" echo in "//trim(echo_file)

 ! To be able to use YAML
 call f_lib_initialize()
 
 !Parsing the input file
 call fdf_init(input_file, echo_file)

 output_file = fdf_string('output', 'sample.out')
 open(newunit=of,file=trim(output_file),action="write")
 call init_random()
 call system%init()
 !--------------TEMP --------------------
 nstates = 1
 nspin = 1
 !---------------------------------------
 call states_init(states, system%basis, nstates, nspin, 1)
 call hamiltonian%init(system, states)
 call scf%init()
 call smear%init()

 close(of)
 call fdf_shutdown() ! no Input after this point


 !SCF loop
 call scf_loop(scf, elsi, hamiltonian, system, states, smear)

 
 !Outputs

 !Release memory
 ! which is not released in final procedure for different types
 call f_lib_finalize()

 !End of the calculation
end program esl_demo
