program esl_demo
  use prec, only : dp,ip
  use iso_fortran_env, only : ou=>OUTPUT_UNIT 
  use system, only : system_t
  use fdf, only : fdf_init, fdf_shutdown, fdf_string
  implicit none

  character(len=100) :: input_file,echo_file,output_file
  type(system_t) :: system
  integer(kind=ip) :: of

  input_file="sample.inp"
  If (command_argument_count() == 1 ) Then
    Call get_command_argument(1, input_file)
  End If
  echo_file=trim(input_file)//".echo"
  write(ou,'(a)')"reading instructions from: "//trim(input_file)//" echo in "//trim(echo_file)
  call fdf_init(input_file, echo_file)
  output_file = fdf_string('output', 'sample.out')
  open(newunit=of,file=trim(output_file),action="write")
  call system%init(of)
  close(of)
  call fdf_shutdown() ! no Input after this point
end program esl_demo
