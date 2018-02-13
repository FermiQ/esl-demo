module esl_ion_interaction_m
  use prec, only : dp
  use esl_constants_m
  use esl_geometry_m

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
  subroutine calculate_periodic(this, geo, forces, eii)
    class(ion_interaction_t), intent(in) :: this
    type(geometry_t),           intent(in) :: geo
    real(dp),           intent(out) :: forces(:,:)
    real(dp),           intent(out) :: eii

    integer :: ia, ja, is, js, ncopy, ic, idim
    integer :: ss, igx, igy, igz, isph
    real(dp) :: dist, zi, zj, rcut, rcopy(3), erfc, charge
    real(dp) :: factor, gx, glen2, gg(3), r(3)
    complex(dp) :: sumat, aa, tmp(3), phase(geo%n_atoms)
    integer :: kk, nbmin(3), nbmax(3)
    eii = 0.d0
    forces(1:3, 1:geo%n_atoms) = 0.d0

    !We have to perform an Ewald summation due to the periodic copies
    rcut = 6.0d0/this%alpha


    !We start with the short range part of the Ewald summation
    do ia = 1, geo%n_atoms
      is = geo%species_idx(ia)
      zi = geo%species(is)%z_ion

      !We need to find the periodic copies with a range of rcut
      gg(1:3) = matmul(geo%xyz(1:3,ia),geo%icell(1:3,1:3))
      do idim = 1, 3
        nbmin(idim) = -nint(-(gg(idim)-rcut)/(geo%cell(idim,idim))+0.5d0)
        nbmax(idim) =  nint((gg(idim)+rcut)/(geo%cell(idim,idim))+0.5d0)
      end do
      ncopy = product(nbmax(1:3)-nbmin(1:3)+1)

      do ic = 1, ncopy
        !get the position of the copy inside rcopy
        kk = ic-1
        do idim = 3, 1, -1
          factor = mod(kk, nbmax(idim)-nbmin(idim)+1) + nbmin(idim)
          rcopy(idim) = gg(idim) - geo%cell(idim,idim)*factor
          kk = kk/(nbmax(idim)-nbmin(idim)+1)
        end do
        rcopy(1:3) = matmul(geo%cell(1:3,1:3), rcopy(1:3))

        do ja = 1, geo%n_atoms
          js = geo%species_idx(ja)
          zj = geo%species(js)%z_ion

          ! Calculate the distance between the two atoms
          r(1:3) = rcopy(1:3) - geo%xyz(1:3,ja)
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
    do ia = 1, geo%n_atoms
      is = geo%species_idx(ia)
      zi = geo%species(is)%z_ion
      charge = charge + zi
      eii = eii - this%alpha/sqrt(PI)*zi**2
    end do

    !Long range part of the Ewald sumation
    rcut = sum(geo%icell(1:3,1)**2)
    do idim = 2, 3
      rcut = min(rcut, sum(geo%icell(1:3,idim)**2))
    end do
    rcut = sqrt(rcut)
    isph = ceiling(9.5d0*this%alpha/rcut)

    !Adding the G=0 term
    eii = eii - PI*charge**2/(2.d0*this%alpha*geo%vol)
    do igx = -isph, isph
      do igy = -isph, isph
        do igz = -isph, isph
          ss = igx**2 + igy**2 + igz**2
          if(ss == 0 .or. ss > isph**2) cycle

          gg(1:3) = igx*geo%icell(1:3,1) + igx*geo%icell(1:3,2) + igx*geo%icell(1:3,3)
          glen2 = sum(gg(1:3)**2)
          !removing te G=0 component
          if (glen2 < 1.0e-10) cycle

          gx = -0.25d0*glen2/this%alpha**2
          if(gx < -36.d0) cycle
          factor = 2*PI/geo%vol*exp(gx)/glen2
          if(factor < epsilon(factor)) cycle

          sumat = cmplx(0.d0,0.d0, kind=dp)
          do ia = 1, geo%n_atoms
            is = geo%species_idx(ia)
            gx = sum(gg(1:3)*geo%xyz(1:3,ia))
            aa = geo%species(is)%z_ion*cmplx(cos(gx),sin(gx),kind=dp)
            phase(ia) = aa
            sumat = sumat + aa
          end do

          eii = eii + real(factor*sumat*conjg(sumat))

          do ia = 1, geo%n_atoms
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
  subroutine calculate_isolated(this, geo, forces, eii)
    class(ion_interaction_t), intent(in)  :: this
    type(geometry_t),         intent(in)  :: geo
    real(dp),                 intent(out) :: forces(:,:)
    real(dp),                 intent(out) :: eii

    integer :: ia, ja, is, js
    real(dp) :: dist, zi, zj, dd, r(3), f(3)

    eii = 0.d0
    forces(1:3, 1:geo%n_atoms) = 0.d0     

    do ia = 1, geo%n_atoms      
      is = geo%species_idx(ia)
      zi = geo%species(is)%z_ion

      do ja = ia + 1, geo%n_atoms 
        js = geo%species_idx(ja)
        zj = geo%species(js)%z_ion

        ! Calculate the distance between the two atoms
        r(1:3) = geo%xyz(1:3,ia) - geo%xyz(1:3,ja)
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
