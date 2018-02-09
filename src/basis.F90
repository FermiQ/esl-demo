module basis_esl

  use prec
  use fdf, only : fdf_get, leqi

  use message_esl

  implicit none
  private

  public ::                &
       basis_t

  type pw_basis_t
     real(dp) :: ecut
     integer :: npw !< Number of plane waves for a given state, spin and kpoint.
  end type pw_basis_t

  !Data structure for the basis
  type basis_t
     integer :: type
     type(pw_basis_t) :: pw_basis
    integer :: size !< Number of coefficients in the basis
    contains
      private
      procedure, public :: init
      procedure, public :: init_basis
      procedure, public :: summary
      final :: cleanup
  end type basis_t

  integer, public, parameter :: &
       PLANEWAVES   = 1, &
       ATOMICORBS   = 2


  contains

   !Initialize the physical system
   !----------------------------------------------------
   subroutine init(this)
     class(basis_t) :: this

     character(len=100) :: str
     real(dp) :: ecut

     str = fdf_get('BasisSet','Planewaves')
     if ( leqi(str, 'Planewaves') )  Then
        this%type = PLANEWAVES
        this%pw_basis%ecut = fdf_get('cut-off', 10._dp, 'Ha')
     else if ( leqi(str, 'AtomicOrbitals') ) then
        this%type = ATOMICORBS
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

   !Initialize the atomic orbitals
   !----------------------------------------------------
   subroutine init_basis(basis, ndims, gcell)
     class(basis_t) :: basis
     integer, dimension(3), intent(in) :: ndims
     real(dp), dimension(3,3), intent(in) :: gcell

     select case (basis%type)
     case (PLANEWAVES)
        basis%pw_basis%npw = getNumberOfPW(ndims, basis%pw_basis%ecut, gcell, [0._dp, 0._dp, 0._dp])
     case (ATOMICORBS)
     end select

   end subroutine init_basis
   
   !Summary
   !----------------------------------------------------
   subroutine summary(basis)
     use yaml_output
     class(basis_t), intent(in) :: basis

     call yaml_mapping_open("basis")
     select case (basis%type)
     case ( PLANEWAVES )
        call yaml_map("Type", "Plane waves")
        call yaml_map("Cut-off (Ha)", basis%pw_basis%ecut)
        call yaml_map("Number of plane-waves", basis%pw_basis%npw)
     case (ATOMICORBS)
        call yaml_map("Type", "Atomic orbitals")
     end select
     call yaml_mapping_close()
   end subroutine summary

   function getNumberOfPW(ndims, ecut, gcell, kpt) result(npw)
     use numerics, only: pi
     implicit none
     integer :: npw
     integer, dimension(3), intent(in) :: ndims
     real(dp), intent(in) :: ecut
     real(dp), dimension(3,3), intent(in) :: gcell
     real(dp), dimension(3), intent(in) :: kpt

     integer :: i1, i2, i3, i
     real(dp), dimension(3,3) :: gmet
     real(dp) :: threshold

     do i = 1, 3
        gmet(i, :) = gcell(1, i) * gcell(1, :) + &
             &   gcell(2, i) * gcell(2, :) + &
             &   gcell(3, i) * gcell(3, :)
     end do

     npw = 0
     threshold = 0.5_dp * ecut / pi**2
     do i1 = -ndims(1) / 2, ndims(1)/2
        do i2 = -ndims(2) / 2, ndims(2)/2
           do i3 = -ndims(3) / 2, ndims(3)/2
              if (dsq(i1, i2, i3) <= threshold) npw = npw + 1
           end do
        end do
     end do
   contains
     function dsq(i1, i2, i3)
       implicit none
       integer, intent(in) :: i1, i2, i3
       real(dp) :: dsq

       dsq=gmet(1,1)*(kpt(1)+dble(i1))**2&
            & +gmet(2,2)*(kpt(2)+dble(i2))**2&
            & +gmet(3,3)*(kpt(3)+dble(i3))**2&
            & +2._dp*(gmet(1,2)*(kpt(1)+dble(i1))*(kpt(2)+dble(i2))&
            & +gmet(2,3)*(kpt(2)+dble(i2))*(kpt(3)+dble(i3))&
            & +gmet(3,1)*(kpt(3)+dble(i3))*(kpt(1)+dble(i1)))

     end function dsq
   end function getNumberOfPW
end module basis_esl
