!< Generecially used precision module
!<
!< All variable kinds should be specified here.
module prec

  implicit none

  public

  !< int precision integer kind
  integer, parameter :: ip = selected_int_kind(9)
  !< Long precision integer kind
  integer, parameter :: lp = selected_int_kind(18)

  !< Single precision real kind
  integer, parameter :: sp = selected_real_kind(6,30)
  !< Double precision real kind
  integer, parameter :: dp = selected_real_kind(14,100)

end module prec
