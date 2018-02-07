module system_esl
  use prec, only : dp,ip
  use numeric, only : invertCell
  use fdf, only : block_fdf, fdf_integer, fdf_block,fdf_defined, &
                   parsed_line, fdf_breals, fdf_bline, fdf_bnames, &
                   fdf_physical

  implicit none
  private

  public ::                &
            system_t,      &
            options_t
  
  !Data structure for the system
  type system_t
    integer(kind=ip) :: nAtoms
    real(kind=dp),allocatable :: coord(:,:) ! (1:3,natoms)
    integer(kind=ip) :: nSpecies
    real(kind=dp) :: cell(3,3)=0.0_dp
    real(kind=dp) :: icell(3,3)=0.0_dp
    character(len=10), dimension(:), allocatable :: el,sp
    character(len=100), dimension(:), allocatable :: potName
    real(dp) :: vol
  contains
    private
    procedure, public :: init
    procedure, public :: summary
    procedure, public :: volume
    final  :: cleanup
  end type system_t

  type :: options_t
    character(len=100) :: output_file
  end type

  contains

   !Initialize the physical system
   !----------------------------------------------------
   subroutine init(sys, of)
     class(system_t),  intent(inout) :: sys
     integer(kind=ip), intent(inout) :: of

     logical :: isdef
     integer :: j,i
     type(block_fdf)            :: blk
     type(parsed_line), pointer :: pline

     isdef = .false.
     sys%nAtoms = fdf_integer('NumberOfAtoms', 0)
     allocate(sys%coord(1:3,sys%nAtoms))
     allocate(sys%el(sys%nAtoms))

     isdef = fdf_defined('coordinates')
     if (isDef) then
       if (fdf_block('coordinates', blk)) then
         j = 1
         do while((fdf_bline(blk, pline)) .and. (j <= sys%nAtoms))
           sys%coord(1:3,j) = [(fdf_breals(pline, i),i=1,3)]
           sys%el(j) = fdf_bnames(pline, 1)
           j = j + 1
         enddo
       endif
     else

     endif

     isDef = fdf_defined('potentials')
     if (isDef) then
       if (fdf_block('potentials', blk)) then
         sys%nSpecies=0
         do while((fdf_bline(blk, pline)))
           sys%nSpecies=sys%nSpecies+1
         enddo
       endif
       allocate(sys%sp(sys%nSpecies),sys%potName(sys%nSpecies))
       if (fdf_block('potentials', blk)) then
         j = 1
         do while((fdf_bline(blk, pline)) .and. (j <= sys%nSpecies))
           sys%sp(j) = fdf_bnames(pline, 1)
           sys%potName(j) = fdf_bnames(pline, 2)
           j = j + 1
         enddo
       endif
     else
     endif

     isdef = fdf_defined('cubic')
     sys%cell=0.0_dp
     if (isDef) then
       sys%cell(1,1)=fdf_physical('cubic', 0.0_dp, 'Ang')
       sys%cell(2,2)=sys%cell(1,1)
       sys%cell(3,3)=sys%cell(1,1)
     endif

     sys%icell=invertCell(sys%cell)
     call sys%summary(of)

   end subroutine init
 
   !Release
   !----------------------------------------------------
   subroutine cleanup(sys)
    class(system_t) :: sys

    if (allocated(sys%coord)) deallocate(sys%coord)
    if (allocated(sys%el)) deallocate(sys%el)
    if (allocated(sys%sp)) deallocate(sys%sp)
    if (allocated(sys%potName)) deallocate(sys%potName)

  end subroutine cleanup

  !Summary
  !----------------------------------------------------
  subroutine summary(sys,of)
    class(system_t) :: sys
    integer(kind=ip) :: of

    integer :: i

    write(of,'(a80)')"Atom Coordinates"
    write(of,'(a11,a20,a20,a20)')"Element |","X|","Y|","Z|"
    do i =1, sys%nAtoms
      write(of,'(a10,3(es19.10,1x))')trim(sys%el(i)),sys%coord(:,i)
    enddo
    write(of,'(a,es16.8,a)')"Volume: ",sys%volume(),"Ang^3"

  end subroutine summary

  !----------------------------------------------------
  real(dp) function volume(sys)
    class(system_t)               :: sys

    volume = sys%cell(1,1)*(sys%cell(2,2)*sys%cell(3,3)-sys%cell(2,3)*sys%cell(3,2)) - &
      sys%cell(1,2)*(sys%cell(2,1)*sys%cell(3,3)-sys%cell(2,3)*sys%cell(3,1)) + &
      sys%cell(1,3)*(sys%cell(2,1)*sys%cell(3,2)-sys%cell(2,2)*sys%cell(3,1))
    sys%Vol = volume

  end function volume


end module system_esl
