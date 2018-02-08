module basis_esl
  use fdf, only : block_fdf, fdf_block,fdf_defined, &
                  parsed_line, fdf_breals, fdf_bline, fdf_bnames, &
                  fdf_get
  use prec

  use message_esl

  implicit none
  private

  public ::                &
       basis_t

  type pw_basis_t
     real(dp) :: ecut
  end type pw_basis_t

  !Data structure for the basis
  type basis_t
     integer :: basis_type
     type(pw_basis_t) :: pw_basis
    integer :: size !< Number of coefficients in the basis
    contains
      private
      procedure, public :: init
      procedure, public :: summary
      final :: cleanup
  end type basis_t

  integer, public, parameter :: &
    PLANEWAVES   = 1,           &
    ATOMICORBS   = 2


  contains

   !Initialize the physical system
   !----------------------------------------------------
   subroutine init(this)
     class(basis_t) :: this

     character(len=100) :: str
     real(dp) :: ecut

     str=fdf_get('BasisSet','Planewaves')
     if (str == 'Planewaves')  Then
        this%basis_type=PLANEWAVES
        this%pw_basis%ecut=fdf_get('cut-off', 10._dp)
     else if (str=='AtomicOrbitals') then
       this%basis_type=ATOMICORBS
     else
       call message_error("The specified basis set is not correct.")
     endif

     !TEMP
     this%size = 1

   end subroutine init
 
   !Release
   !----------------------------------------------------
   subroutine cleanup(this)
     type(basis_t) :: this

   end subroutine cleanup

   subroutine summary(basis)
     use yaml_output
     implicit none
     class(basis_t), intent(in) :: basis

     call yaml_mapping_open("basis")
     select case (basis%basis_type)
     case (PLANEWAVES)
        call yaml_map("Type", "Plane waves")
        call yaml_map("Cut-off (Ha)", basis%pw_basis%ecut)
     case (ATOMICORBS)
        call yaml_map("Type", "Atomic orbitals")
     end select
     call yaml_mapping_close()
   end subroutine summary

end module basis_esl
