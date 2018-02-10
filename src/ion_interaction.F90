module esl_ion_interaction_m

  use prec, only : dp

  use esl_system_m

  implicit none
  private

  public :: ion_interaction_t

  !Data structure for the ion interaction
  type ion_interaction_t
    real(dp) :: alpha !< Ewald sum parameter
  contains

    procedure, public :: init
    procedure, public :: calculate_periodic
    procedure, public :: calculate_isolated
    final  :: cleanup

  end type ion_interaction_t

contains

  !Initialize 
  !----------------------------------------------------
  subroutine init(this)
    class(ion_interaction_t), intent(inout) :: this

  end subroutine init

  !Release
  !----------------------------------------------------
  subroutine cleanup(this)
    type(ion_interaction_t), intent(inout) :: this

  end subroutine cleanup

  !Calc the ion-ion interaction for periodic systems
  !----------------------------------------------------
  subroutine calculate_periodic(this, sys, forces, eii)
    use esl_constants_m, only : PI

    class(ion_interaction_t), intent(in) :: this
    type(system_t),           intent(in) :: sys
    real(dp),           intent(out) :: forces(:,:)
    real(dp),           intent(out) :: eii

    integer :: ia, ja, is, js, ncopy, ic, idim
    integer :: ss, igx, igy, igz, isph
    real(dp) :: dist, zi, zj, rcut, rcopy(3), erfc, charge
    real(dp) :: factor, gx, glen2, gg(3), r(3)
    complex(dp) :: sumat, aa, tmp(3), phase(sys%geo%n_atoms)

    eii = 0.d0
    forces(1:3, 1:sys%geo%n_atoms) = 0.d0

    !We have to perform an Ewald summation due to the periodic copies
    rcut = 6.0d0/this%alpha


    !We start with the short range part of the Ewald summation
    do ia = 1, sys%geo%n_atoms
      is = sys%geo%species_idx(ia)
      zi = sys%geo%species(is)%Z

      !We need to find the periodic copies with a range of rcut
      ncopy = 1
      do ic = 1, ncopy
        !get the position of the copy inside rcopy

        do ja = 1, sys%geo%n_atoms
          js = sys%geo%species_idx(ja)
          zj = sys%geo%species(js)%Z

          ! Calculate the distance between the two atoms
          r(1:3) = rcopy(1:3) - sys%geo%xyz(1:3,ja)
          dist = sum(r(1:3) ** 2) ** 0.5_dp

          erfc = 1.d0 - erf(this%alpha*dist)

          !the force
          forces(1:3,ia) = forces(1:3,ia) &
              - zi*zj*r(1:3)*(erfc/dist + 2.d0*this%alpha/sqrt(PI)*exp(-(this%alpha*dist)**2))/dist**2

          !energy
          eii = eii + 0.5d0*zi*zj*erfc/dist
        end do
      end do
    end do

    !Self interaction part
    do ia = 1, sys%geo%n_atoms
      is = sys%geo%species_idx(ia)
      zi = sys%geo%species(is)%Z
      charge = charge + zi
      eii = eii - this%alpha/sqrt(PI)*zi**2
    end do

    !Long range part of the Ewald sumation
    rcut = sum(sys%geo%icell(1:3,1)**2)
    do idim = 2, 3
      rcut = min(rcut, sum(sys%geo%icell(1:3,idim)**2))
    end do
    rcut = sqrt(rcut)
    isph = ceiling(9.5d0*this%alpha/rcut)

    !Adding the G=0 term
    eii = eii - PI*charge**2/(2.d0*this%alpha*sys%geo%vol)
    do igx = -isph, isph
      do igy = -isph, isph
        do igz = -isph, isph
          ss = igx**2 + igy**2 + igz**2
          if(ss == 0 .or. ss > isph**2) cycle

          gg(1:3) = igx*sys%geo%icell(1:3,1) + igx*sys%geo%icell(1:3,2) + igx*sys%geo%icell(1:3,3)
          glen2 = sum(gg(1:3)**2)
          !removing te G=0 component
          if(glen2 < 1.0e-10) cycle

          gx = -0.25d0*glen2/this%alpha**2
          if(gx < -36.d0) cycle
          factor = 2*PI/sys%geo%vol*exp(gx)/glen2
          if(factor < epsilon(factor)) cycle

          sumat = cmplx(0.d0,0.d0, kind=dp)
          do ia = 1, sys%geo%n_atoms
            is = sys%geo%species_idx(ia)
            gx = sum(gg(1:3)*sys%geo%xyz(1:3,ia))
            aa = sys%geo%species(is)%Z*cmplx(cos(gx),sin(gx))
            phase(ia) = aa
            sumat = sumat + aa
          end do

          eii = eii + real(factor*sumat*conjg(sumat))

          do ia = 1, sys%geo%n_atoms
            tmp(1:3) = cmplx(0.d0,1.d0,kind=dp)*gg(1:3)*phase(ia)
            forces(1:3,ia) = forces(1:3,ia) &
                -factor*real(conjg(tmp)*sumat+tmp*conjg(tmp))
          end do

        end do
      end do
    end do

  end subroutine calculate_periodic

  !Calc the ion-ion interaction for isolated systems
  !----------------------------------------------------
  subroutine calculate_isolated(this, sys, forces, eii)
    class(ion_interaction_t), intent(in) :: this
    type(system_t),           intent(in) :: sys
    real(dp),           intent(out) :: forces(:,:)
    real(dp),           intent(out) :: eii

    integer :: ia, ja, is, js
    real(dp) :: dist, zi, zj, dd, r(3), f(3)

    eii = 0.d0
    forces(1:3, 1:sys%geo%n_atoms) = 0.d0     

    do ia = 1, sys%geo%n_atoms      
      is = sys%geo%species_idx(ia)
      zi = sys%geo%species(is)%Z

      do ja = ia + 1, sys%geo%n_atoms 
        js = sys%geo%species_idx(ja)
        zj = sys%geo%species(js)%Z

        ! Calculate the distance between the two atoms
        r(1:3) = sys%geo%xyz(1:3,ia) - sys%geo%xyz(1:3,ja)
        dist = sum(r(1:3) ** 2) ** 0.5_dp

        !the force
        dd = zi*zj/dist
        f(1:3) = (dd/dist**2)*r(1:3)
        forces(1:3,ia) = forces(1:3,ia) + f(1:3)
        forces(1:3,ja) = forces(1:3,ja) - f(1:3)

        !energy
        eii = eii + dd

      end do
    end do

  end subroutine calculate_isolated

end module esl_ion_interaction_m
