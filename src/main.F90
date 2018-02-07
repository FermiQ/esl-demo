program esl_demo
 use basis_esl
 use hamiltonian_esl
 use parser_esl
 use scf_esl
 use system_esl

 implicit none

 type(basis_t)       :: basis
 type(hamiltonian_t) :: hamiltonian
 type(parser_t)      :: parser
 type(system_t)      :: system
 type(scf_t)         :: scf

 !Init MPI

 !Init data basis strucutes 

 !Read Input variables and init corresponding data structures
 call parser_init(parser)

 call system_init(system, parser)

 call basis_init(basis)

 call hamiltonian_init(hamiltonian, parser)

 call scf_init(scf, parser)

 !SCF loop
 call scf_loop(scf)

 
 !Outputs

 !Release memory
 call scf_end(scf)

 call hamiltonian_end(hamiltonian)

 call basis_end(basis)

 call system_end(system)

 call parser_end(parser) 

 !End of the calculation
 
end program
