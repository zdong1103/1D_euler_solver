module M_parameter

  implicit none

  integer, parameter :: IP  = selected_int_kind(8)
  integer, parameter :: SP  = 4
  integer, parameter :: DP  = 8
  real,    parameter :: inf = huge(DP)

  type, public :: material_information
    real(DP)     :: gamma
    real(DP)     :: rho_s
    real(DP)     :: R  = 8.3144598_DP
    real(DP)     :: PI = 3.141592654_DP
  end type material_information

end module M_parameter