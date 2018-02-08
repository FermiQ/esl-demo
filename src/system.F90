module system_esl
  use prec, only : dp,ip
  use numeric_esl, only : invert_cell
  use fdf, only : block_fdf, fdf_integer, fdf_block,fdf_defined, &
                   parsed_line, fdf_breals, fdf_bline, fdf_bnames, &
                   fdf_physical

  use basis_esl
  use grid_esl
  use smear_esl
  use states_esl

  implicit none
  private

  public ::                &
            system_t
  
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

    type(basis_t) :: basis
    type(grid_t)  :: grid
    real(dp) :: nElectrons
  contains
    private
    procedure, public :: init
    procedure, public :: summary
    procedure, public :: volume
    final  :: cleanup
  end type system_t

  contains

   !Initialize the physical system
   !----------------------------------------------------
    subroutine init(sys)
     class(system_t) :: sys

     logical :: isdef
     integer :: j,i
     type(block_fdf)            :: blk
     type(parsed_line), pointer :: pline

     integer :: nstates, nspin

     call sys%basis%init()

     isdef = fdf_defined('cubic')
     sys%cell=0.0_dp
     if (isDef) then
       sys%cell(1,1)=fdf_physical('cubic', 0.0_dp, 'Bohr')
       sys%cell(2,2)=sys%cell(1,1)
       sys%cell(3,3)=sys%cell(1,1)
     endif
     sys%icell=invert_cell(sys%cell)

     !Init the grid
     call sys%grid%init(sys%basis, sys%cell(1,1))

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

     call sys%summary()

   end subroutine init
 
   !Release
   !----------------------------------------------------
   subroutine cleanup(sys)
     type(system_t) :: sys

    if (allocated(sys%coord)) deallocate(sys%coord)
    if (allocated(sys%el)) deallocate(sys%el)
    if (allocated(sys%sp)) deallocate(sys%sp)
    if (allocated(sys%potName)) deallocate(sys%potName)

  end subroutine cleanup

  !Summary
  !----------------------------------------------------
  subroutine summary(sys)
    use yaml_output
    class(system_t) :: sys

    integer :: i

    call yaml_mapping_open("System")
    call yaml_map("Cell", sys%cell)
    call yaml_sequence_open("Atom Coordinates", advance = "no")
    call yaml_comment("Element | X| Y| Z|", hfill = "-")
    do i =1, sys%nAtoms
       call yaml_sequence(advance="no")
       call yaml_map(trim(sys%el(i)), sys%coord(:,i))
    enddo
    call yaml_sequence_close()
    call yaml_map("Volume (Bohr^3)", sys%volume())
    call yaml_mapping_close()

    call sys%grid%summary()
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
