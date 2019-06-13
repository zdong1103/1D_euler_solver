module M_integration
  ! ======================= Include Modules =======================
  use M_precision
  use M_pin
  use M_domain
  use M_spatial_recon

  ! ========================= Declarations ========================
  implicit none

  ! ============================= Body ============================
  contains
    !! ---------------------- calculate dt ------------------------

    subroutine calc_dt(myd, pin, dt)
      implicit none
      type(parameter_input), intent(in), target  :: pin
      type(domain), intent(in)  :: myd  ! my_domain
      real(DP),     intent(out) :: dt
      real(DP)      :: a = 0._DP, a_new ! sound speed
      integer(IP)   :: i

      do i = myd%is, myd%ie
        call sound_speed(myd%cells(i)%cons(1), myd%cells(i)%cons(2), myd%cells(i)%cons(3), &
                        &myd%material_info%gamma, a_new)
        if ( a_new .gt. a) then
          a = a_new
        end if
      end do

      dt = pin%CFL * myd%L / myd%N / a

    end subroutine calc_dt

    !! ----------------- integrate dt using RK2 -------------------

    subroutine integrate_step(myd, dt)
      type(domain), intent(inout)       :: myd
      type(domain)                      :: pseudo
      type(flux)                        :: F
      real(DP), dimension(:,:), pointer :: k1
      real(DP), dimension(:,:), pointer :: k2
      real(DP), intent(in)    :: dt
      integer(IP) :: i, j

      allocate(k1(3,myd%N))
      allocate(k2(3,myd%N))

      ! k1 = dt * f(u,t)
      call calc_flux(F, myd)
      call copy_domain(pseudo, myd)

      variable_loop_k1: do j = 1, 3
        domain_loop_k1: do i = 1, myd%N
          k1(j,i) = dt * (- myd%L / myd%N ) * (F%F_R(j,i) - F%F_L(j,i))
          pseudo%cells(i+myd%nghost)%cons(j) = pseudo%cells(i+myd%nghost)%cons(j) &
                                           & + k1(j,i) / 2._DP
          call cons2prim(pseudo%cells(i+myd%nghost), myd%material_info%gamma)
        end do domain_loop_k1
      end do variable_loop_k1

      ! update ghost cells
      do i = 1, myd%nghost
        call copy_cell_info(pseudo%cells(i), pseudo%cells(1+myd%nghost))
        call copy_cell_info(pseudo%cells(i+myd%N), pseudo%cells(myd%N+myd%nghost))
      end do

      call clear_flux(F)

      ! k2 = dt * f(u+k1/2,t+dt/2)
      call calc_flux(F, pseudo)

      variable_loop_k2: do j = 1, 3
        domain_loop_k2: do i = 1, myd%N
          k2(j,i) = dt * (- myd%L / myd%N ) * (F%F_R(j,i) - F%F_L(j,i))
          myd%cells(i+myd%nghost)%cons(j) = myd%cells(i+myd%nghost)%cons(j) &
                                        & + k2(j,i)
          call cons2prim(myd%cells(i+myd%nghost), myd%material_info%gamma)
        end do domain_loop_k2
      end do variable_loop_k2

      ! update ghost cells
      do i = 1, myd%nghost
        call copy_cell_info(myd%cells(i), myd%cells(1+myd%nghost))
        call copy_cell_info(myd%cells(i+myd%N+myd%nghost), myd%cells(myd%N+myd%nghost))
      end do

      myd%t = myd%t + dt
      myd%nstep = myd%nstep + 1_IP

      deallocate(k1)
      deallocate(k2)
      call clear_flux(F)
      call delete_domain(pseudo)

    end subroutine integrate_step

end module M_integration

