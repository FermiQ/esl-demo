module message_esl

  implicit none
  private

  public ::                   &
       message_error

contains

  !Print an error message and leave the code
  !----------------------------------------------------
  subroutine message_error(str)
    character(len=*)  :: str

    print *, str
    stop

  end subroutine message_error

end module message_esl
