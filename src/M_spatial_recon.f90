module M_spatial_recon
  ! ======================= Include Modules =======================
  use M_precision
  use M_domain

  ! ========================= Declarations ========================
  implicit none
  
  type, public :: flux_param
    ! For all conservative variables
    real(DP), dimension(:,:), pointer :: phi     ! phi(r_i) is the flux limiter
    real(DP), dimension(:,:), pointer :: u_RR    ! u_(i+1/2)^R
    real(DP), dimension(:,:), pointer :: u_RL    ! u_(i+1/2)^L
    real(DP), dimension(:,:), pointer :: u_LR    ! u_(i-1/2)^R
    real(DP), dimension(:,:), pointer :: u_LL    ! u_(i-1/2)^L

  end type flux_param

  type, public :: flux
    ! For all conservative variables
    real(DP), dimension(:),   pointer :: a_R   ! a_(i+1/2)
    real(DP), dimension(:),   pointer :: a_L   ! a_(i-1/2)   
    real(DP), dimension(:,:), pointer :: F_R   ! F_(i+1/2)*
    real(DP), dimension(:,:), pointer :: F_L   ! F_(i-1/2)*

  end type flux
  ! ============================= Body ============================
  contains
    !! --------- initialize and calculate flux parameters ---------

    subroutine calc_flux_param(this, myd)
      implicit none
      type(flux_param), intent(inout) :: this
      type(domain),     intent(in)    :: myd   ! my_domain
      integer(IP) :: i, j
      real(DP)    :: ri                        ! r_i

      allocate(this%phi(3,myd%N_tot))
      allocate(this%u_RR(3,myd%N_tot))
      allocate(this%u_RL(3,myd%N_tot))
      allocate(this%u_LR(3,myd%N_tot))
      allocate(this%u_LL(3,myd%N_tot))

      ! ------------------------------------------------------------
      ! r_i = (u_i - u_i-1) / (u_i+1 - u_i)
      ! phi(r_i) = max[0, min(1, r_i)]
      ! ------------------------------------------------------------
      variable_loop_1 : do j = 1, 3
        domain_loop_1 : do i = myd%is, myd%ie
          if ( myd%cells(i+1)%cons(j) .eq. myd%cells(i)%cons(j) ) then
            if ( myd%cells(i)%cons(j) .ge. myd%cells(i-1)%cons(j) ) then
              ri = inf
            else
              ri = -1 * inf
            endif
          else
            ri = (myd%cells(i)%cons(j) - myd%cells(i-1)%cons(j)) / &
                &(myd%cells(i+1)%cons(j) - myd%cells(i)%cons(j))
          end if

          ! minmod flux limiter
          this%phi(j,i) = max(0._DP, min(1._DP, ri))
          write (*,*) "phi (", j, ",", i, ") ", this%phi(j,i)
        end do domain_loop_1
      end do variable_loop_1
      
      ! ------------------------------------------------------------
      ! u_(i+1/2)^R = u_i+1 - 0.5 * phi(r_i+1) * (u_i+2 - u_i+1)
      ! u_(i+1/2)^L = u_i   + 0.5 * phi(r_i)   * (u_i+1 - u_i)
      ! u_(i-1/2)^R = u_i   - 0.5 * phi(r_i)   * (u_i+1 - u_i)
      ! u_(i-1/2)^L = u_i-1 + 0.5 * phi(r_i-1) * (u_i   - u_i-1)
      ! ------------------------------------------------------------
      variable_loop_2 : do j = 1, 3
        domain_loop_2 : do i = myd%is, myd%ie
          this%u_RR(j,i) = myd%cells(i+1)%cons(j) - 0.5_DP * this%phi(j,i+1) * &
                         &(myd%cells(i+2)%cons(j) - myd%cells(i+1)%cons(j))
          
          this%u_RL(j,i) = myd%cells(i)%cons(j)   + 0.5_DP * this%phi(j,i)   * &
                         &(myd%cells(i+1)%cons(j) - myd%cells(i)%cons(j))
          
          this%u_LR(j,i) = myd%cells(i)%cons(j)   - 0.5_DP * this%phi(j,i)   * &
                         &(myd%cells(i+1)%cons(j) - myd%cells(i)%cons(j))
          
          this%u_LL(j,i) = myd%cells(i-1)%cons(j) + 0.5_DP * this%phi(j,i)   * &
                         &(myd%cells(i)%cons(j)   - myd%cells(i-1)%cons(j))          
        end do domain_loop_2
      end do variable_loop_2

    end subroutine calc_flux_param

    !! ----------------------- Destroyer --------------------------

    subroutine clear_flux_param(this)
      implicit none
      type(flux_param), intent(inout) :: this
      
      deallocate(this%phi)
      deallocate(this%u_RR)
      deallocate(this%u_RL)
      deallocate(this%u_LR)
      deallocate(this%u_LL)

    end subroutine clear_flux_param

    !! ----------------- calculate momentum flux ------------------

    subroutine momentum_flux(rho, rhoU, E, gamma, MomFlux)
      implicit none
      real(DP), intent(in)  :: rho, rhoU, E, gamma
      real(DP), intent(out) :: MomFlux
      
      ! conservation of momentum
      MomFlux = 0.5_DP * (3._DP - gamma) * rhoU * rhoU / rho + (gamma - 1._DP) * E
    end subroutine momentum_flux

    !! ------------------ calculate energy flux -------------------

    subroutine energy_flux(rho, rhoU, E, gamma, EFlux)
      implicit none
      real(DP), intent(in)  :: rho, rhoU, E, gamma
      real(DP), intent(out) :: EFlux

      !conservation of energy
      EFlux = gamma * rhoU * E / rho - 0.5_DP * (gamma - 1._DP) * rhoU * rhoU * rhoU / rho / rho
    end subroutine energy_flux

    !! ------------------ calculate sound speed -------------------

    subroutine sound_speed(rho, rhoU, E, gamma, a)
      implicit none
      real(DP), intent(in)  :: rho, rhoU, E, gamma
      real(DP), intent(out) :: a

      a = sqrt(gamma * (gamma - 1._DP) * (E - 0.5_DP * rhoU * rhoU / rho) / rho)
    end subroutine sound_speed

    !! -------------------- calculate fluxes ----------------------

    subroutine calc_flux(this, myd)
      type(flux),       intent(inout) :: this
      type(domain),     intent(in)    :: myd    ! my_domain
      type(flux_param)                :: myf    ! my_flux_param
      integer(IP) :: i, j
      real(DP)    :: gamma
      real(DP)    :: a_RR, a_RL, a_LR, a_LL
      real(DP)    :: F_RR, F_RL, F_LR, F_LL     ! F(u_(i+1/2)^R), F(u_(i+1/2)^L)
                                                ! F(u_(i-1/2)^R), F(u_(i-1/2)^L)
      
      gamma = myd%material_info%gamma
      call calc_flux_param(myf,myd)

      allocate(this%a_R(myd%N_tot))
      allocate(this%a_L(myd%N_tot))
      allocate(this%F_R(3,myd%N_tot))
      allocate(this%F_L(3,myd%N_tot))
      write (*,*) "-----------------------------------------------"

      ! ------------------------------------------------------------
      ! a_(i+-1/2) = max[eigenvalue(dF_(i+-1/2)^R/du), eigenvalue(dF_(i+-1/2)^L/du)]
      ! eigenvalues for F are U, U-a, U+a, the max one is U+a
      ! a (sound speed) = sqrt(gamma*p/rho)
      ! ------------------------------------------------------------
      do i = myd%is, myd%ie
        call sound_speed(myf%u_RR(1,i), myf%u_RR(2,i), myf%u_RR(3,i), gamma, a_RR)
        call sound_speed(myf%u_RL(1,i), myf%u_RL(2,i), myf%u_RL(3,i), gamma, a_RL)
        call sound_speed(myf%u_LR(1,i), myf%u_LR(2,i), myf%u_LR(3,i), gamma, a_LR)
        call sound_speed(myf%u_LL(1,i), myf%u_LL(2,i), myf%u_LL(3,i), gamma, a_LL)

        this%a_R(i) = max(myf%u_RR(2,i)/myf%u_RR(1,i) + a_RR, myf%u_RL(2,i)/myf%u_RL(1,i) + a_RL)
        this%a_L(i) = max(myf%u_LR(2,i)/myf%u_LR(1,i) + a_LR, myf%u_LL(2,i)/myf%u_LL(1,i) + a_LL)
        write (*,*) "a_R, a_L (", i, ") ", this%a_R(i), this%a_L(i)
      end do

      variable_loop : do j = 1, 3
        domain_loop : do i = myd%is, myd%ie
          select case (j) ! (1,2,3) => (rho, rhoU, E)
          case (1)
            F_RR = myf%u_RR(1,i)
            F_RL = myf%u_RL(1,i)
            F_LR = myf%u_LR(1,i)
            F_LL = myf%u_LL(1,i)
          case (2)
            call momentum_flux(myf%u_RR(1,i), myf%u_RR(2,i), myf%u_RR(3,i), gamma, F_RR)
            call momentum_flux(myf%u_RL(1,i), myf%u_RL(2,i), myf%u_RL(3,i), gamma, F_RL)
            call momentum_flux(myf%u_LR(1,i), myf%u_LR(2,i), myf%u_LR(3,i), gamma, F_LR)
            call momentum_flux(myf%u_LL(1,i), myf%u_LL(2,i), myf%u_LL(3,i), gamma, F_LL)
          case (3)
            call energy_flux(myf%u_RR(1,i), myf%u_RR(2,i), myf%u_RR(3,i), gamma, F_RR)
            call energy_flux(myf%u_RL(1,i), myf%u_RL(2,i), myf%u_RL(3,i), gamma, F_RL)
            call energy_flux(myf%u_LR(1,i), myf%u_LR(2,i), myf%u_LR(3,i), gamma, F_LR)
            call energy_flux(myf%u_LL(1,i), myf%u_LL(2,i), myf%u_LL(3,i), gamma, F_LL)
          end select

          ! ------------------------------------------------------------
          ! F_(1+1/2)* = 0.5 * (F(u_(i+1/2)^R) + F(u_(i+1/2)^L)) 
          !              - a_(i+1/2) * (u_(i+1/2)^R - u_(i+1/2)^L)
          ! F_(1-1/2)* = 0.5 * (F(u_(i-1/2)^R) + F(u_(i-1/2)^L)) 
          !              - a_(i-1/2) * (u_(i-1/2)^R - u_(i-1/2)^L)
          ! ------------------------------------------------------------
          this%F_R(j,i) = 0.5_DP * ((F_RR + F_RL) - this%a_R(i) * (myf%u_RR(j,i) - myf%u_RL(j,i)))
          this%F_L(j,i) = 0.5_DP * ((F_LR + F_LL) - this%a_L(i) * (myf%u_LR(j,i) - myf%u_LL(j,i)))
          ! write (*,*) "F_R, F_L (", j, ",", i, ") ", this%F_R(j,i), this%F_L(j,i)
        end do domain_loop
      end do variable_loop

      call clear_flux_param(myf)

    end subroutine calc_flux

    !! ----------------------- Destroyer --------------------------

    subroutine clear_flux(this)
      implicit none
      type(flux), intent(inout) :: this

      deallocate(this%a_R)
      deallocate(this%a_L)
      deallocate(this%F_R)
      deallocate(this%F_L)

    end subroutine clear_flux

end module M_spatial_recon

