program main
  ! ======================= Include Modules =======================
  use M_precision
  use M_domain
  use M_pin
  
  ! ========================= Declarations ========================
  implicit none

  character(len=30) :: filename = "input.txt"
  type(parameter_input), target :: pin
  type(domain) :: my_domain

  ! ============================= Body ============================
  call init_read(pin, filename)
  call init_domain(my_domain, pin)
  
  call delete_domain(my_domain)

  write (*,*) "success!"

end program  