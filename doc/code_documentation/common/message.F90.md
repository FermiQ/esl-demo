## Overview

The `esl_message_m` module offers a straightforward utility for handling fatal errors within the application. It provides a single public subroutine, `message_error`, which is designed to print a user-specified error message to the standard output and then immediately terminate the program's execution. This serves as a simple, global error-handling mechanism.

## Key Components

- **Module:** `esl_message_m`
    - **Description:** A utility module focused on providing a standardized way to report a fatal error message and stop program execution.
- **Public Subroutine:**
    - `message_error(str)`: The sole public procedure in this module. It prints the provided error message and halts the program.

### Subroutine/Function: `message_error(str)`
- **Description:** This subroutine handles fatal error conditions. When called:
    1.  It prints the input character string `str` to the standard output device (typically the console, using `print *`).
    2.  It then immediately terminates the entire program's execution using the Fortran `stop` statement.
    This provides a clear indication to the user that an unrecoverable error has occurred and why.
- **Arguments:**
    - `str (character(len=*))`: The error message string that will be displayed to the user just before the program halts. The `len=*` indicates that it can accept strings of any length.
- **Returns:** Not applicable (Subroutine, as the program terminates within it).

## Important Variables/Constants
- This module does not define any public Fortran `parameter` constants.

## Usage Examples
```fortran
! TODO: Add usage example
! Conceptual usage within another part of the code:

module data_processor_m
  use esl_message_m, only: message_error
  implicit none

  public :: process_data

contains
  subroutine process_data(input_value)
    integer, intent(in) :: input_value

    if (input_value < 0) then
      call message_error("Error: Input value cannot be negative in process_data.")
      ! Program execution stops here if the condition is met.
    else if (input_value == 0) then
      call message_error("Fatal: Zero is not a permitted value for processing.")
    end if

    ! ... (proceed with normal data processing if no error) ...
    print *, "Processing data with value: ", input_value

  end subroutine process_data
end module data_processor_m

! program test_data_processor
!   use data_processor_m
!   call process_data(10)  ! OK
!   call process_data(-5)  ! Will call message_error and stop.
! end program test_data_processor
```

## Dependencies and Interactions

- **Internal Dependencies:**
    - This module relies only on standard Fortran intrinsic procedures (`print *`, `stop`). It does not have dependencies on other custom modules within the ESL library for its core functionality.
- **External Libraries:**
    - This module has no dependencies on external libraries.
- **Interactions with other components:**
    - **Error Handling Strategy:** `message_error` is intended to be invoked by any other module or subroutine within the application when a critical, non-recoverable error is detected. Instead of implementing complex error propagation mechanisms (like returning error codes up the call stack), a call to `message_error` provides a direct and immediate way to halt execution.
    - **User Notification:** Its primary interaction is to provide a final message to the user (or to a log file, if standard output is redirected) explaining the reason for the program's termination.
    - **Debugging:** In development, calls to `message_error` can help pinpoint the location and cause of fatal errors quickly.
    - **Simplicity:** It offers a very simple API for a common task, reducing boilerplate error-handling code in other parts of the application for unrecoverable situations.
