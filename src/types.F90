module types 
  use prec, only : dp,ip
  implicit none
  private

  type, public :: systemT
    integer(kind=ip) :: nAtoms
    real(kind=dp),allocatable :: x(:),y(:),z(:)
    integer(kind=ip) :: nSpecies
    real(kind=dp) :: cell(3,3)=0.0_dp
    real(kind=dp) :: icell(3,3)=0.0_dp
    character(len=10), dimension(:), allocatable :: el,sp
    character(len=100), dimension(:), allocatable :: potName
    real(dp) :: vol
  contains 
    private
    procedure, public :: summary
    procedure, public :: volume
    final  :: cleanup

  end type

  type, public :: optionsT
    character(len=100) :: output_file
  end type
contains
  subroutine cleanup(sys)
    type(systemT) :: sys
    if (allocated(sys%x)) deallocate(sys%x) 
    if (allocated(sys%y)) deallocate(sys%y) 
    if (allocated(sys%z)) deallocate(sys%z)
    if (allocated(sys%el)) deallocate(sys%el)
    if (allocated(sys%sp)) deallocate(sys%sp)
    if (allocated(sys%potName)) deallocate(sys%potName)

  end subroutine cleanup
  subroutine summary(sys,of)
  class(systemT) :: sys
    integer(kind=ip) :: of

    integer :: i

    write(of,'(a80)')"Atom Coordinates"
    write(of,'(a11,a20,a20,a20)')"Element |","X|","Y|","Z|"
    do i =1, sys%nAtoms
      write(of,'(a10,3(es19.10,1x))')trim(sys%el(i)),sys%x(i),sys%y(i),sys%z(i)
    enddo 
    write(of,'(a,es16.8,a)')"Volume: ",sys%volume(),"Ang^3"
  end subroutine summary
  real(dp) function volume(sys)
  class(systemT)               :: sys

    volume = sys%cell(1,1)*(sys%cell(2,2)*sys%cell(3,3)-sys%cell(2,3)*sys%cell(3,2)) - &
      sys%cell(1,2)*(sys%cell(2,1)*sys%cell(3,3)-sys%cell(2,3)*sys%cell(3,1)) + &
      sys%cell(1,3)*(sys%cell(2,1)*sys%cell(3,2)-sys%cell(2,2)*sys%cell(3,1)) 
    sys%Vol = volume
  end function volume
end module types
