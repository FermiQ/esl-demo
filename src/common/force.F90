!< Force module implementing a type to retain the different force contributions
!<
!< The currently collected forces are:
!<
!<  - local contributions
!<  - non-local contribitions
!<  - ion-ion interactions
!<
!< All forces are in Hartree/Bohr
module esl_force_m

  use prec, only : dp

  implicit none
  private

  public :: force_t

  !< Force data-structure to retain forces
  type force_t

    !< Sum of all force-terms (total force)
    real(dp), allocatable :: total(:,:)
    !< Local contributions
    real(dp), allocatable :: loc(:,:)
    !< Non-local contributions
    real(dp), allocatable :: nl(:,:)
    !< Ion-ion interactions
    real(dp), allocatable :: ionion(:,:)

  contains
    private

    procedure, public :: init
    procedure, public :: calculate
    final :: cleanup

  end type force_t

contains

  !< Initialize the force object with `natoms`
  subroutine init(this, natoms)
    class(force_t), intent(inout) :: this
    integer, intent(in) :: natoms

    allocate(this%total(3, natoms))
    this%total(:, :) = 0._dp
    allocate(this%loc(3, natoms))
    this%loc(:, :) = 0._dp
    allocate(this%nl(3, natoms))
    this%nl(:, :) = 0._dp
    allocate(this%ionion(3, natoms))
    this%ionion(:, :) = 0._dp

  end subroutine init

  !< Calculate and update the total forces
  subroutine calculate(this)
    class(force_t), intent(inout) :: this
    this%total(:, :) = this%loc(:, :) + this%nl(:, :) + this%ionion(:, :)
  end subroutine calculate

  !< Print to std-out the total force
  !<
  !< @param[in] only_total if true, only print-out the total force, else everything is written.
  subroutine display(this, only_total)
    use yaml_output
    class(force_t) :: this
    logical, intent(in), optional :: only_total
    
    logical :: lonly_total

    lonly_total = .true.
    if ( present(only_total) ) &
        lonly_total = only_total

    if ( lonly_total ) then
      call yaml_map("Force total", this%total)
    else
      call yaml_mapping_open("Forces", advance='NO')
      call yaml_comment("Hartree/Bohr", hfill = "-")
      call yaml_map("Total", this%total)
      call yaml_map("Local", this%loc)
      call yaml_map("Non Local", this%nl)
      call yaml_map("Ion-Ion", this%ionion)
    end if
    call yaml_mapping_close()

  end subroutine display
  
  !< Clean up the object
  subroutine cleanup(this)
    type(force_t), intent(inout) :: this

    if ( allocated(this%total) ) &
        deallocate(this%total)
    if ( allocated(this%loc) ) &
        deallocate(this%loc)
    if ( allocated(this%nl) ) &
        deallocate(this%nl)
    if ( allocated(this%ionion) ) &
        deallocate(this%ionion)

  end subroutine cleanup

end module esl_force_m
