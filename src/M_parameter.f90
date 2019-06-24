! Written by Zhenyang DONG

module M_parameter

! ========================= Declarations ========================
  implicit none

  integer,  parameter :: IP  = selected_int_kind(8)
  integer,  parameter :: SP  = 4
  integer,  parameter :: DP  = 8
  real,     parameter :: inf = huge(DP)
  real(DP), parameter :: R  = 8.31446261815324
  real(DP), parameter :: PI = 3.14159265358979

  type, public :: material_information
    real(DP)     :: gamma
    real(DP)     :: rho_s
  end type material_information

  ! ============================= Body ============================
  contains
    !! ------------------ calculate sound speed -------------------

    subroutine sound_speed(rho, rhoU, E, mif, c)
      implicit none
      type(material_information), intent(in)  :: mif
      real(DP), intent(in)  :: rho, rhoU, E
      real(DP), intent(out) :: c

      c = sqrt(mif%gamma * (mif%gamma - 1._DP) * (E - 0.5_DP * rhoU * rhoU / rho) / rho)
    end subroutine sound_speed

end module M_parameter