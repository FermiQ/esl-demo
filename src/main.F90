program esl_demo
 use scf_esl

 implicit none

 !Init MPI

 !Init data basis strucutes 

 !Read Input variables and init corresponding data structures

 !SCF
 call scf_init()

 call scf_loop()

 call scf_end()

 !Outputs

 !Release memory

 !End of the calculation

 
end program
