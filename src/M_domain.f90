module M_domain
  ! ======================= Include Modules =======================
  use M_precision
  use M_pin

  ! ======================== Declarations =========================
  implicit none
  
  type, public :: cell
    ! Information about in this cell
    integer(IP)  :: index
    real(DP)     :: x

    ! conservative variables
    real(DP)     :: rho       ! density
    real(DP)     :: rhoU      ! momentum
    real(DP)     :: E         ! energy

    ! primitive variables
    real(DP)     :: p         ! pressure
    real(DP)     :: U         ! velocity

  end type cell
  
  type, public :: domain
    ! Information about the entire domain
    real(DP)     :: t         ! time
    real(DP)     :: L         ! domain length
    integer(IP)  :: N         ! total # cells
    integer(IP)  :: nstep     ! step #

    type(material_information), pointer  :: material_info   ! general gas parameters
    type(cell), dimension(:),   pointer  :: cells           ! cell array

  end type domain

  ! ======================== Subroutines =========================
  contains
    subroutine init_cell(this, index, pin) ! Initialize cell information
      implicit none
      type(cell),            intent(inout)       :: this
      integer(IP),           intent(in)          :: index
      type(parameter_input), intent(in), target  :: pin

      this%index = index
      this%x     = (index - 0.5_DP) * pin%L / pin%N
      if ( index .le. pin%N/2._DP ) then
        this%rho = pin%rho_L
        this%p   = pin%p_L
        this%U   = pin%U_L
      else
        this%rho = pin%rho_R
        this%p   = pin%p_R
        this%U   = pin%U_R
      end if
      this%rhoU = this%rho * this%U
      this%E    = this%p / (pin%material_info%gamma - 1._DP) + 0.5_DP * this%rho * this%U * this%U
    end subroutine init_cell

    subroutine delete_cell(this)
      implicit none
      type(cell), intent(inout) :: this

      this%index = 0_IP
      this%rho   = 0._DP
      this%p     = 0._DP
      this%U     = 0._DP
      this%rhoU  = 0._DP
      this%E     = 0._DP
    end subroutine delete_cell

    subroutine init_domain(this, pin)
      implicit none
      type(domain),          intent(inout)       :: this
      type(parameter_input), intent(in), target  :: pin
      integer(IP) :: i

      this%N = pin%N
      this%L = pin%L
      this%t = 0._DP
      this%nstep = 0_IP
      this%material_info => pin%material_info
      allocate(this%cells(pin%N))

      do i = 1, pin%N
        call init_cell(this%cells(i),i, pin)
      end do
    end subroutine init_domain

    subroutine delete_domain(this)
      implicit none
      type(domain), intent(inout) :: this
      integer(IP) :: i

      do i = 1, this%N
        call delete_cell(this%cells(i))
      end do

      deallocate(this%cells)
      this%N = 0_IP
      this%L = 0._DP
      this%t = 0._DP
      this%nstep = 0_IP
    end subroutine delete_domain

end module M_domain
