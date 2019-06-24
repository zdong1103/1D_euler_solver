! Written by Zhenyang DONG

program main
  ! ======================= Include Modules =======================
  use M_parameter
  use M_domain
  use M_pin
  use M_output
  use M_spatial_recon
  use M_integration
  
  ! ========================= Declarations ========================
  implicit none

  character(len=30) :: filename = "input.txt"
  type(parameter_input), target :: pin
  real(DP)     :: Time  ! total time
  real(DP)     :: dt    ! Time step
  type(domain) :: my_domain

  ! ============================= Body ============================
  call init_read(pin, filename)
  call init_domain(my_domain, pin)
  call output(my_domain)
  Time = pin%t_tot

  do while ( my_domain%t .lt. Time )
    call calc_dt(my_domain, pin, dt)
    ! adjust dt of the last step
    if ( my_domain%t + dt .gt. Time ) then
      dt = Time - my_domain%t
    end if

    write(*,"(A, I3, A, E9.3, A, E9.3)") "cycle = ", my_domain%nstep+1, "  dt = ", dt, "  t = ", my_domain%t
    call integrate_step(my_domain, dt)
    call output(my_domain)
  end do

  call delete_domain(my_domain)

  write (*,*) "Complete!"

end program  
