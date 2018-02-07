program esl_demo
  use prec, only : dp,ip
  use iso_fortran_env, only : ou=>OUTPUT_UNIT 
  use types, only : systemT,optionsT
  use parser, only : parse_input
  use useful, only : invertCell
  implicit none

  character(len=100) :: input_file
  type(optionsT) :: options
  type(systemT) :: system
  integer(kind=ip) :: of

  input_file="sample.inp"
  If (command_argument_count() == 1 ) Then
    Call get_command_argument(1, input_file)
  End If
  write(ou,'(a)')"reading instructions from: "//trim(input_file)
  call parse_input(input_file,system,options)
  open(newunit=of,file=trim(options%output_file),action="write")
  system%icell=invertCell(system%cell)
  call system%summary(of)
  close(of)
end program esl_demo
