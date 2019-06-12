module M_precision

  implicit none

  integer, parameter :: IP  = selected_int_kind(8)
  integer, parameter :: SP  = 4
  integer, parameter :: DP  = 8
  real,    parameter :: inf = huge(DP)

end module M_precision