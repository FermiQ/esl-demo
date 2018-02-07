module parser
  use prec
  use types, only : systemT, optionsT
  USE fdf
  implicit none
  private

  public :: parse_input
contains
  subroutine parse_input(filename,system,options)
    character(len=*), intent(in) :: filename
    type(systemT), intent(out) :: system
    type(optionsT), intent(out) :: options

    character(len=100) :: echofile
    logical :: isdef
    integer :: j
    type(block_fdf)            :: blk
    type(parsed_line), pointer :: pline


    isdef = .false.
    echofile=trim(filename)//".echo"
    call fdf_init(filename, echofile)
    system%nAtoms = fdf_integer('NumberOfAtoms', 0)
    allocate(system%x(system%nAtoms),system%y(system%nAtoms),system%z(system%nAtoms))
    allocate(system%el(system%nAtoms))

    options%output_file = fdf_string('output', 'sample.out')
    isdef = fdf_defined('coordinates')
    if (isDef) then
      if (fdf_block('coordinates', blk)) then
        j = 1
        do while((fdf_bline(blk, pline)) .and. (j <= system%nAtoms))
          system%x(j) = fdf_breals(pline, 1)
          system%y(j) = fdf_breals(pline, 2)
          system%z(j) = fdf_breals(pline, 3)
          system%el(j) = fdf_bnames(pline, 1)
          j = j + 1
        enddo
      endif

    else
    endif
    isDef = fdf_defined('potentials')
    if (isDef) then
      if (fdf_block('potentials', blk)) then
        system%nSpecies=0
        do while((fdf_bline(blk, pline)))
          system%nSpecies=system%nSpecies+1
        enddo

      endif
      allocate(system%sp(system%nSpecies),system%potName(system%nSpecies))
      if (fdf_block('potentials', blk)) then
        j = 1
        do while((fdf_bline(blk, pline)) .and. (j <= system%nSpecies))
          system%sp(j) = fdf_bnames(pline, 1)
          system%potName(j) = fdf_bnames(pline, 2)
          j = j + 1
        enddo
      endif

    else
    endif
    isdef = fdf_defined('cubic')
    system%cell=0.0_dp
    if (isDef) then
      system%cell(1,1)=fdf_physical('cubic', 0.0_dp, 'Ang')
      system%cell(2,2)=system%cell(1,1)
      system%cell(3,3)=system%cell(1,1)
    endif

    call fdf_shutdown()
  end subroutine parse_input
end module parser

