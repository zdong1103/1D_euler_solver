module M_output
  ! ======================= Include Modules =======================
  use M_precision
  use M_domain

  ! ======================== Declarations =========================
  implicit none

  character(len=20) :: fmt_var = "(2A8,5A12)"
  character(len=20) :: fmt_num = "(2F8.3,5F12.3)"

  ! ======================== Subroutines =========================
  contains
    subroutine output(this)
      implicit none
      type(domain), intent(inout) :: this
      character(len=20) :: output_file
      character(len=5)  :: frame_num
      integer(IP)       :: i
      logical           :: exist

      if ( this%nstep .le. 10 ) then
        write(frame_num,"(A,I1)") '0000',this%nstep
      elseif ( this%nstep .le. 100 ) then
        write(frame_num,"(A,I2)") '000',this%nstep
      elseif ( this%nstep .le. 1000 ) then
        write(frame_num,"(A,I3)") '00',this%nstep
      elseif ( this%nstep .le. 10000 ) then
        write(frame_num,"(A,I4)") '0',this%nstep
      else
        write(frame_num,"(I5)") this%nstep
      end if
      
      output_file = "output/out."//frame_num//".dat"
      inquire(file=output_file, exist=exist)

      if ( exist ) then
        open(1, file=output_file, status='old',action="write")
      else
        open(1, file=output_file, status='new',action="write")
      end if
      
      write(1, fmt_var) "t","x","rho","U","p","rhoU","E"
      do i = 1, this%N
        write(1,fmt_num) this%t, this%cells(i)%x, this%cells(i)%rho, this%cells(i)%U, &
                         &this%cells(i)%p, this%cells(i)%rhoU, this%cells(i)%E
      end do
      close(1)

    end subroutine output

end module M_output