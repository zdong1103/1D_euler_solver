module M_pin
  ! ======================= Include Modules =======================
  use M_precision

  ! ======================== Declarations =========================
  implicit none

  type, public :: material_information
    real(DP)     :: gamma
    real(DP)     :: R  = 8.3144598_DP
    real(DP)     :: PI = 3.141592654_DP
  end type material_information

  type, public :: parameter_input
    ! Domain size
    integer(IP)  :: N
    integer(IP)  :: nghost    ! # ghost cells on one side
    real(DP)     :: L

    ! Numerical parameter
    real(DP)     :: t_tot
    real(DP)     :: CFL

    ! Initial conditions
    ! Left side
    real(DP)     :: p_L         ! pressure
    real(DP)     :: U_L         ! velocity
    real(DP)     :: rho_L       ! density

    ! Right side
    real(DP)     :: p_R         ! pressure
    real(DP)     :: U_R         ! velocity
    real(DP)     :: rho_R       ! density

    type(material_information) :: material_info

    ! contains
    !   procedure  :: init   => init_read

  end type parameter_input

  ! ======================== Subroutines =========================
  contains
    subroutine init_read(this, filename)
      implicit none
      type(parameter_input), intent(inout) :: this
      character(len=30),     intent(in)    :: filename
      character(len=30),     dimension(20) :: lines
      character(len=30)  :: var_name
      character(len=30)  :: value
      integer(IP) :: istat, i, index

      open(1, file=filename, status='old')
      i = 1
      do while(.true.)
        read(1, "(A)", iostat=istat) lines(i)
        lines(i) = trim(lines(i))
        index     = scan(lines(i),"=")
        var_name  = trim(lines(i)(1:index-1))
        value     = trim(lines(i)(index+1:))
        select case (var_name)
        case ("N") 
          read(value,*) this%N
        case ("nghost")
          read(value,*) this%nghost
        case ("L") 
          read(value,*) this%L
        case ("t") 
          read(value,*) this%t_tot
        case ("CFL")
          read(value,*) this%CFL
        case ("p_L")
          read(value,*) this%p_L
        case ("U_L")
          read(value,*) this%U_L
        case ("rho_L")
          read(value,*) this%rho_L
        case ("p_R")
          read(value,*) this%p_R
        case ("U_R")
          read(value,*) this%U_R
        case ("rho_R")
          read(value,*) this%rho_R
        case ("gamma")
          read(value,*) this%material_info%gamma
        end select
        if ( istat /= 0 ) exit
        i = i + 1
      end do
      close(1)
    end subroutine init_read

end module M_pin

