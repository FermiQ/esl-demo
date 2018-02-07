program esl_demo
 use basis_esl, only : basis_T
 use fdf, only : fdf_init, fdf_shutdown, fdf_string
 use hamiltonian_esl, only : hamiltonian_t
 use iso_fortran_env, only : ou=>OUTPUT_UNIT
 use prec, only : dp, ip
 use scf_esl, only : scf_t, scf_loop
 use system_esl, only : system_t
 use numeric_esl, only : init_random

 implicit none

 type(basis_t)       :: basis
 type(hamiltonian_t) :: hamiltonian
 type(system_t)      :: system
 type(scf_t)         :: scf
 character(len=100) :: input_file,echo_file,output_file
 integer :: of
 !Init MPI

 !Init data basis strucutes 

 !Read Input variables and init corresponding data structures
 input_file="sample.inp"
 if (command_argument_count() == 1 ) Then
   call get_command_argument(1, input_file)
 end if

 echo_file=trim(input_file)//".echo"
 write(ou,'(a)')"reading instructions from: "//trim(input_file)//" echo in "//trim(echo_file)
 
 !Parsing the input file
 call fdf_init(input_file, echo_file)

 output_file = fdf_string('output', 'sample.out')
 open(newunit=of,file=trim(output_file),action="write")
 call init_random()
 call system%init(of)
 call basis%init()
 call hamiltonian%init(system)
 call scf%init()

 close(of)
 call fdf_shutdown() ! no Input after this point


 !SCF loop
 call scf_loop(scf)

 
 !Outputs

 !Release memory
 ! which is not released in final procedure for different types

 !End of the calculation
end program esl_demo
