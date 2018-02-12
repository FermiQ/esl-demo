module esl_geometry_m
  use prec
  use fdf

  use esl_numeric_m, only : matr3inv
  use esl_species_m
  use esl_message_m

  implicit none
  private

  public ::                &
       geometry_t

  !Data structure for the system
  type geometry_t
    ! Species
    integer(ip) :: n_species
    type(species_t), allocatable :: species(:)

    ! Atoms
    integer(ip) :: n_atoms
    real(dp),    allocatable :: xyz(:,:) ! (1:3,1:natoms)
    integer(ip), allocatable :: species_idx(:) ! (1:n_atoms)
    
    ! Cell
    real(dp) :: cell(3,3) = 0.0_dp
    real(dp) :: icell(3,3) = 0.0_dp
    real(dp) :: vol

  contains
    private
    procedure, public :: init
    procedure, public :: summary
    procedure, public :: volume
    final  :: cleanup
  end type geometry_t

contains

  !Initialize the geometry
  !----------------------------------------------------
  subroutine init(this)
    class(geometry_t) :: this

    logical :: is_def
    integer :: is, ia, i
    type(block_fdf)            :: blk
    type(parsed_line), pointer :: pline
    
    ! Read cell input options
    is_def = fdf_defined('cubic')
    this%cell = 0.0_dp
    if (is_def) then
      this%cell(1,1) = fdf_get('cubic', 0.0_dp, 'Bohr')
      this%cell(2,2) = this%cell(1,1)
      this%cell(3,3) = this%cell(1,1)
    endif
    this%icell = matr3inv(this%cell)

    ! Read species options
    this%n_species = 0
    is_def = fdf_defined('species')
    if (is_def) then
      if (fdf_block('species', blk)) then
        do while ((fdf_bline(blk, pline)))
          this%n_species = this%n_species + 1
        end do
      endif
      allocate(this%species(this%n_species))
      if (fdf_block('species', blk)) then
        is = 1
        do while((fdf_bline(blk, pline)) .and. (is <= this%n_species))
          call this%species(is)%init(trim(fdf_bnames(pline, 1)), fdf_bnames(pline, 2))
          is = is + 1
        end do
      end if
    endif

    if (this%n_species == 0) then
      call message_error("No species defined in input file!")    
    end if

    ! Read atoms options
    is_def = .false.
    this%n_atoms = fdf_get('NumberOfAtoms', 0)
    allocate(this%xyz(1:3, this%n_atoms))
    allocate(this%species_idx(this%n_atoms))

    is_def = fdf_defined('coordinates')
    if (is_def) then
      if (fdf_block('coordinates', blk)) then
        ia = 1
        do while((fdf_bline(blk, pline)) .and. (ia <= this%n_atoms))
          this%xyz(1:3, ia) = [(fdf_breals(pline, i), i=1,3)]

          do is = 1, this%n_species + 1
            if (is == this%n_species + 1) &
              call message_error("Species "//trim(fdf_bnames(pline, 1))//" is unknown!")
            if (leqi(fdf_bnames(pline, 1), this%species(is)%label)) then
              this%species_idx(ia) = is
              exit
            end if
          end do
          ia = ia + 1            
        end do
      end if
    else
      call message_error("No atomic coordinates defined in input file!")
    end if

    call this%summary()

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(geometry_t) :: this

    integer :: i

    if (allocated(this%species)) deallocate(this%species)
    if (allocated(this%xyz)) deallocate(this%xyz)
    if (allocated(this%species_idx)) deallocate(this%species_idx)

  end subroutine cleanup

  !Summary
  !----------------------------------------------------
  subroutine summary(this)
    use yaml_output
    class(geometry_t) :: this

    integer :: ia

    call yaml_mapping_open("Geometry")
    call yaml_map("Cell", this%cell)
    call yaml_sequence_open("Atom Coordinates", advance = "no")
    call yaml_comment("Element | X| Y| Z|", hfill = "-")
    do ia = 1, this%n_atoms
       call yaml_sequence(advance="no")
       call yaml_map(trim(this%species(this%species_idx(ia))%label), this%xyz(:,ia))
    enddo
    call yaml_sequence_close()
    call yaml_map("Volume (Bohr^3)", this%volume())
    call yaml_mapping_close()

  end subroutine summary

  !----------------------------------------------------
  function volume(this) result(vol)
    class(geometry_t), intent(inout) :: this
    real(dp) :: vol

    vol = this%cell(1,1)*(this%cell(2,2)*this%cell(3,3)-this%cell(2,3)*this%cell(3,2)) - &
          this%cell(1,2)*(this%cell(2,1)*this%cell(3,3)-this%cell(2,3)*this%cell(3,1)) + &
          this%cell(1,3)*(this%cell(2,1)*this%cell(3,2)-this%cell(2,2)*this%cell(3,1))
    this%Vol = vol

  end function volume

end module esl_geometry_m
