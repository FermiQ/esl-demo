module esl_basis_pw_m

  use prec
  use yaml_output
  use esl_geometry_m
  
  private

  public :: basis_pw_t

  type basis_pw_t
    real(dp) :: ecut !< Plane wave cut-off in Hartree
    integer  :: npw  !< Number of plane waves
  contains
    private
    procedure, public :: init
    procedure, public :: summary
    final :: cleanup
  end type basis_pw_t
  
contains

  !Initialize the basis
  !----------------------------------------------------
  subroutine init(this, ecut, ndims, gcell)
    class(basis_pw_t) :: this
    real(dp), intent(in) :: ecut
    integer,  intent(in) :: ndims(3)
    real(dp), intent(in) :: gcell(3,3)

    this%ecut = ecut
    this%npw = get_number_of_pw(ndims, ecut, gcell, [0._dp, 0._dp, 0._dp])

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(basis_pw_t) :: this
  end subroutine cleanup

  !Summary
  !----------------------------------------------------
  subroutine summary(this)
    class(basis_pw_t), intent(in) :: this

    call yaml_mapping_open("basis_pw")
    call yaml_map("Cut-off (Ha)", this%ecut)
    call yaml_map("Number of plane-waves", this%npw)
    call yaml_mapping_close()

  end subroutine summary

  integer function get_number_of_pw(ndims, ecut, gcell, kpt) result(npw)
    use esl_constants_m, only: pi
    integer,  intent(in) :: ndims(3)
    real(dp), intent(in) :: ecut
    real(dp), intent(in) :: gcell(3, 3)
    real(dp), intent(in) :: kpt(3)

    integer :: i1, i2, i3, i
    real(dp) :: gmet(3, 3)
    real(dp) :: threshold

    do i = 1, 3
       gmet(i, :) = gcell(1, i) * gcell(1, :) &
            &        + gcell(2, i) * gcell(2, :) &
            &        + gcell(3, i) * gcell(3, :)
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
      integer, intent(in) :: i1, i2, i3
      real(dp) :: dsq

      dsq = gmet(1, 1)*(kpt(1) + dble(i1))**2 &
           & + gmet(2, 2)*(kpt(2) + dble(i2))**2 &
           & + gmet(3, 3)*(kpt(3) + dble(i3))**2 &
           & + 2._dp*(gmet(1, 2)*(kpt(1) + dble(i1))*(kpt(2) + dble(i2)) &
           & + gmet(2, 3)*(kpt(2) + dble(i2))*(kpt(3) + dble(i3)) &
           & + gmet(3, 1)*(kpt(3) + dble(i3))*(kpt(1) + dble(i1)))

    end function dsq
    
  end function get_number_of_pw

end module esl_basis_pw_m
