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
    real(DP), dimension(3) :: cons  ! (1,2,3) => (rho, rhoU, E)

    ! primitive variables
    real(DP), dimension(2) :: prim  ! (1,2) => (p, U)

  end type cell
  
  type, public :: domain
    ! Information about the entire domain
    real(DP)     :: t         ! time
    real(DP)     :: L         ! domain length
    integer(IP)  :: N         ! # cells
    integer(IP)  :: nghost    ! # ghost cells on each side
    integer(IP)  :: is        ! starting cell index
    integer(IP)  :: ie        ! ending cell index
    integer(IP)  :: nstep     ! step #

    type(material_information), pointer  :: material_info   ! general gas parameters
    type(cell), dimension(:),   pointer  :: cells           ! cell array

  end type domain

  ! ======================== Subroutines =========================
  contains
    !! ------------------- initialize a cell ---------------------

    subroutine init_cell(this, index, pin) ! Initialize cell information
      implicit none
      type(cell),            intent(inout)       :: this
      integer(IP),           intent(in)          :: index
      type(parameter_input), intent(in), target  :: pin

      this%index = index
      this%x     = (index - 0.5_DP - pin%nghost) * pin%L / pin%N
      if ( this%x .le. pin%L/2._DP ) then
        this%cons(1) = pin%rho_L
        this%prim(1) = pin%p_L
        this%prim(2) = pin%U_L
      else
        this%cons(1) = pin%rho_R
        this%prim(1) = pin%p_R
        this%prim(2) = pin%U_R
      end if
      this%cons(2) = this%cons(1) * this%prim(2)
      this%cons(3) = this%prim(1) / (pin%material_info%gamma - 1._DP) + 0.5_DP * &
                    &this%cons(1) * this%prim(2) * this%prim(2)
    end subroutine init_cell

    !! ---------------------- copy a cell ------------------------

    subroutine copy_cell(this, anothercell)
      type(cell),  intent(inout) :: this
      type(cell),  intent(in)    :: anothercell

      this%cons(1) = anothercell%cons(1)
      this%cons(2) = anothercell%cons(2)
      this%cons(3) = anothercell%cons(3)
      this%prim(1) = anothercell%prim(1)
      this%prim(2) = anothercell%prim(2)      
    end subroutine copy_cell

    !! ----------------------- destroyer -------------------------

    subroutine delete_cell(this)
      implicit none
      type(cell), intent(inout) :: this

      this%index   = 0_IP
      this%x       = 0._DP
      this%cons(1) = 0._DP
      this%cons(2) = 0._DP
      this%cons(3) = 0._DP
      this%prim(1) = 0._DP
      this%prim(2) = 0._DP
    end subroutine delete_cell

    !! ------------------- initialize domain ---------------------

    subroutine init_domain(this, pin)
      implicit none
      type(domain),          intent(inout)       :: this
      type(parameter_input), intent(in), target  :: pin
      integer(IP) :: i
      integer(IP) :: N_tot

      this%N      = pin%N
      this%nghost = pin%nghost
      this%is     = pin%nghost + 1_IP
      this%ie     = this%is + this%N - 1_IP
      this%L      = pin%L
      this%t      = 0._DP
      this%nstep  = 0_IP
      this%material_info => pin%material_info

      N_tot = pin%N + 2_IP * pin%nghost
      allocate(this%cells(N_tot))

      do i = 1, N_tot
        call init_cell(this%cells(i),i, pin)
      end do
    end subroutine init_domain

    !! ----------------------- destroyer -------------------------

    subroutine delete_domain(this)
      implicit none
      type(domain), intent(inout) :: this
      integer(IP) :: i
      integer(IP) :: N_tot

      N_tot = this%N + 2_IP * this%nghost
      do i = 1, N_tot
        call delete_cell(this%cells(i))
      end do

      deallocate(this%cells)
      this%N      = 0_IP
      this%nghost = 0_IP
      this%is     = 0_IP
      this%ie     = 0_IP
      this%L      = 0._DP
      this%t      = 0._DP
      this%nstep  = 0_IP
    end subroutine delete_domain

end module M_domain
