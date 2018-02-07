module smear_esl

  implicit none
  private

  public ::                          &
            smear_t,                 &
            smear_calc_fermi,        &
            smear_calc_occ
  
  !Data structure for the system
  type smear_t
    integer :: smearing !< which smearing do we use
  end type smear_t

   integer, public, parameter :: &
    COLD     = 1,           &
    FD       = 2,           &
    MP       = 3,           &
    GAUSSIAN = 4


  contains

   !Compute the Fermi level
   !----------------------------------------------------
   subroutine smear_calc_fermi(this)
     type(smear_t), intent(inout) :: this

     !TODO: Use ELSI here

     select case(this%smearing)
       case(COLD)
       case(FD)
       case(MP)
       case(GAUSSIAN)
     end select

   end subroutine smear_calc_fermi
 
   !Compute the occupations
   !----------------------------------------------------
   subroutine smear_calc_occ(this)
     type(smear_t), intent(inout) :: this

     !TODO: Use ELSI here

     select case(this%smearing)
       case(COLD)
       case(FD)
       case(MP)
       case(GAUSSIAN)
     end select

   end subroutine smear_calc_occ

end module smear_esl
