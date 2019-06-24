! Written by Zhenyang DONG

module M_domain
  ! ======================= Include Modules =======================
  use M_parameter
  use M_pin

  ! ======================== Declarations =========================
  implicit none
  
  type, public :: cell
    ! Information about in this cell
    integer(IP)  :: index
    real(DP)     :: x

    ! conservative variables
    ! (1,2,3) => (rho, rhoU, E)
    real(DP), dimension(3) :: cons  

    ! primitive variables
    ! (1,2) => (p,U)
    real(DP), dimension(2) :: prim  

  end type cell
  
  type, public :: domain
    ! Information about the entire domain
    real(DP)     :: t         ! time
    real(DP)     :: L         ! domain length
    integer(IP)  :: N         ! # cells
    integer(IP)  :: nghost    ! # ghost cells on each side
    integer(IP)  :: N_tot     ! # cells + 2 * # ghost cells
    integer(IP)  :: is        ! starting cell index
    integer(IP)  :: ie        ! ending cell index
    integer(IP)  :: nstep     ! step #

    type(material_information)           :: mif   ! general gas parameters
    type(cell), dimension(:),   pointer  :: cells ! cell array

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
      this%cons(3) = this%prim(1) / (pin%mif%gamma - 1._DP) + 0.5_DP * &
                    &this%cons(1) * this%prim(2) * this%prim(2)
    end subroutine init_cell

    !! ---------------------- copy a cell ------------------------

    subroutine copy_cell(this, anothercell)
      implicit none
      type(cell),  intent(inout) :: this
      type(cell),  intent(in)    :: anothercell

      this%index   = anothercell%index
      this%x       = anothercell%x
      this%cons(1) = anothercell%cons(1)
      this%cons(2) = anothercell%cons(2)
      this%cons(3) = anothercell%cons(3)
      this%prim(1) = anothercell%prim(1)
      this%prim(2) = anothercell%prim(2)      
    end subroutine copy_cell

    subroutine copy_cell_info(this, anothercell)
      implicit none
      type(cell),  intent(inout) :: this
      type(cell),  intent(in)    :: anothercell

      this%cons(1) = anothercell%cons(1)
      this%cons(2) = anothercell%cons(2)
      this%cons(3) = anothercell%cons(3)
      this%prim(1) = anothercell%prim(1)
      this%prim(2) = anothercell%prim(2)
    end subroutine copy_cell_info

    !! --------------- conservative to primitive -----------------

    subroutine cons2prim(this, mif)
      implicit none
      type(cell), intent(inout) :: this
      type(material_information), intent(in)  :: mif

      this%prim(1) = (mif%gamma - 1._DP) * (this%cons(3) - 0.5_DP * this%cons(2) * &
                   &  this%cons(2) / this%cons(1))
      this%prim(2) = this%cons(2) / this%cons(1)

    end subroutine cons2prim

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
      this%N_tot  = pin%N + 2_IP * pin%nghost
      this%is     = pin%nghost + 1_IP
      this%ie     = this%is + this%N - 1_IP
      this%L      = pin%L
      this%t      = 0._DP
      this%nstep  = 0_IP
      this%mif%gamma = pin%mif%gamma

      allocate(this%cells(this%N_tot))

      do i = 1, this%N_tot
        call init_cell(this%cells(i),i, pin)
      end do
    end subroutine init_domain

    !! ------------------- create domain copy --------------------

    subroutine copy_domain(this, anotherdomain)
      implicit none
      type(domain), intent(in)    :: anotherdomain
      type(domain), intent(inout) :: this
      integer(IP) :: i

      this%N      = anotherdomain%N
      this%nghost = anotherdomain%nghost
      this%N_tot  = anotherdomain%N_tot
      this%is     = anotherdomain%is
      this%ie     = anotherdomain%ie
      this%L      = anotherdomain%L
      this%t      = anotherdomain%t
      this%nstep  = anotherdomain%nstep
      this%mif%gamma = anotherdomain%mif%gamma

      allocate(this%cells(this%N_tot))

      do i = 1, this%N_tot
        call copy_cell(this%cells(i), anotherdomain%cells(i))
      end do

    end subroutine copy_domain

    !! ----------------------- destroyer -------------------------

    subroutine delete_domain(this)
      implicit none
      type(domain), intent(inout) :: this
      integer(IP) :: i

      do i = 1, this%N_tot
        call delete_cell(this%cells(i))
      end do

      deallocate(this%cells)
      this%N      = 0_IP
      this%nghost = 0_IP
      this%N_tot  = 0_IP
      this%is     = 0_IP
      this%ie     = 0_IP
      this%L      = 0._DP
      this%t      = 0._DP
      this%nstep  = 0_IP
      this%mif%gamma = 0._DP

    end subroutine delete_domain

end module M_domain
