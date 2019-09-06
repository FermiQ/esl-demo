!< Module implementing a geometry container
!<
!< This geometry is the basic container for atoms and their
!< atomic species. The unit cell is also contained in this
!< type.
!<
!< All lengths are in Bohr.
module esl_geometry_m
  use prec
  use fdf

  use esl_numeric_m, only : matr3inv
  use esl_species_m
  use esl_message_m

  implicit none
  private

  public :: geometry_t

  !< Geometry data structure
  !<
  !< This data structure retains information such as:
  !<  - Number of atoms
  !<  - Atomic coordinates
  !<  - Atomic species (and their specific species type)
  !<  - Cell vectors
  !<  - Cell volume
  type geometry_t
    !< Number of different species
    integer(ip) :: n_species = 0
    !< All different species used in this geometry
    type(species_t), allocatable :: species(:)

    !< Number of atoms
    integer(ip) :: n_atoms = 0
    !< Atomic coordinates, size `3, n_atoms`
    real(dp), allocatable :: xyz(:,:)
    !< Atomic specie indices (referring to `species`), size `n_atoms`
    integer(ip), allocatable :: species_idx(:)
    
    !< Unit cell
    real(dp) :: cell(3,3) = 0.0_dp
    !< Inverse unit cell TODO are the inverse vectors corresponding to the same indicies in cell, or the transposed?
    real(dp) :: icell(3,3) = 0.0_dp
    !< Unit cell volume
    real(dp) :: vol = 0._dp

  contains
    private
    
    procedure, public :: init
    procedure, public :: summary
    procedure, public :: volume
    procedure, public :: electronic_charge
    final  :: finalizer
    
  end type geometry_t

contains

  !< Initialize the geometry by reading input from the user
  subroutine init(this)
    class(geometry_t) :: this

    integer :: is, ia, i
    type(block_fdf)            :: blk
    type(parsed_line), pointer :: pline
    logical :: cell_defined
    
    ! Read cell input options
    cell_defined = fdf_defined('Cubic')
    if ( cell_defined ) then
      
      this%cell(1,1) = fdf_get('Cubic', 0.0_dp, 'Bohr')
      this%cell(2,2) = this%cell(1,1)
      this%cell(3,3) = this%cell(1,1)
      
    end if

    ! Read species options
    this%n_species = 0
    if ( fdf_block('Species', blk) ) then
      do while ( fdf_bline(blk, pline) )
        this%n_species = this%n_species + 1
      end do
    end if
    
    allocate(this%species(this%n_species))
    if ( fdf_block('Species', blk) ) then
      is = 1
      do while ( fdf_bline(blk, pline) .and. is <= this%n_species )
        call this%species(is)%init(trim(fdf_bnames(pline, 1)), fdf_bnames(pline, 2))
        is = is + 1
      end do
    end if

    if ( this%n_species == 0 ) then
      call message_error("No species defined in input file!")    
    end if

    ! Read atoms options
    ! We don't have to put numberofatoms, so we simply count lines
    this%n_atoms = 0
    if ( fdf_block('Coordinates', blk) ) then
      do while ( fdf_bline(blk, pline) )
        this%n_atoms = this%n_atoms + 1
      end do
    end if

    ! Get user-defined number of atoms
    this%n_atoms = fdf_get('NumberOfAtoms', this%n_atoms)
    
    allocate(this%xyz(3, this%n_atoms))
    allocate(this%species_idx(this%n_atoms))
    if ( fdf_block('Coordinates', blk) ) then
      ia = 1
      do while ( fdf_bline(blk, pline) .and. ia <= this%n_atoms )
        this%xyz(1, ia) = fdf_breals(pline, 1)
        this%xyz(2, ia) = fdf_breals(pline, 2)
        this%xyz(3, ia) = fdf_breals(pline, 3)

        do is = 1, this%n_species + 1
          if ( is == this%n_species + 1 ) &
              call message_error("Species "//trim(fdf_bnames(pline, 1))//" is unknown!")
          if ( leqi(fdf_bnames(pline, 1), this%species(is)%label) ) then
            this%species_idx(ia) = is
            exit
          end if
        end do
        ia = ia + 1            
      end do
    else
      call message_error("No atomic coordinates defined in input file!")
    end if

    if ( .not. cell_defined ) then
      
      ! Make it a molecule
      this%cell = 0._dp
      this%cell(1,1) = maxval(this%xyz(1, :), dim=1) - minval(this%xyz(1, :), dim=1)
      this%cell(2,2) = maxval(this%xyz(2, :), dim=1) - minval(this%xyz(2, :), dim=1)
      this%cell(3,3) = maxval(this%xyz(3, :), dim=1) - minval(this%xyz(3, :), dim=1)

      ! Add vacuum of ~ 30 Bohr
      do i = 1, 3
        this%cell(i,i) = this%cell(i,i) + 30._dp
      end do
      
    end if

    ! Ensure that the volume is calculated and stored in the geometry
    if ( this%volume() < 0._dp ) then
      call message_error("Cell volume is negative!")
    end if
    this%icell = matr3inv(this%cell)

  end subroutine init

  !< Clean up the object
  subroutine finalizer(this)
    type(geometry_t) :: this

    if ( allocated(this%species) ) &
        deallocate(this%species)
    if ( allocated(this%xyz) ) &
        deallocate(this%xyz)
    if ( allocated(this%species_idx) ) &
        deallocate(this%species_idx)

  end subroutine finalizer

  !< Print out the geometry coordinates and species information in the YAML output
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
    end do
    call yaml_sequence_close()
    call yaml_sequence_open("Species info")
    do ia = 1, this%n_species
      call yaml_sequence(advance="no")
      call this%species(ia)%summary()
    end do
    call yaml_sequence_close()
    call yaml_map("Number of electrons", this%electronic_charge())
    call yaml_map("Volume (Bohr^3)", this%volume())
    call yaml_mapping_close()

  end subroutine summary

  !< Calculate the volume of the unit-cell
  function volume(this) result(vol)
    class(geometry_t), intent(inout) :: this
    real(dp) :: vol

    vol = this%cell(1,1)*(this%cell(2,2)*this%cell(3,3)-this%cell(2,3)*this%cell(3,2)) - &
          this%cell(1,2)*(this%cell(2,1)*this%cell(3,3)-this%cell(2,3)*this%cell(3,1)) + &
          this%cell(1,3)*(this%cell(2,1)*this%cell(3,2)-this%cell(2,2)*this%cell(3,1))

    ! Since some lattice vectors *may* be "negative" we need to take the absolute value
    ! It rarely happens though.
    this%Vol = abs(vol)

  end function volume

  !< Return the valence electrons for the atoms in this geometry
  function electronic_charge(this) result(charge)
    class(geometry_t), intent(inout) :: this
    real(dp) :: charge

    integer :: ia

    charge = 0.0_dp
    do ia = 1, this%n_atoms
      charge = charge + this%species(this%species_idx(ia))%q
    end do

  end function electronic_charge
  
end module esl_geometry_m
