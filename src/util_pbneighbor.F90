! util_pbneighbor
!> @brief Calc the neighbor vector in periodic boundary condition

subroutine util_pbneighbor(vx, imirror)

  use const_maxsize
  use var_inp, only : inperi

  implicit none

  ! --------------------------------------------------------------------
  real(PREC), intent(inout) :: vx(3)
  integer, intent(out) :: imirror

  ! --------------------------------------------------------------------
  ! local variables
  integer :: ix
  integer :: imi(3)

  ! --------------------------------------------------------------------

  do ix = 1, 3
     if(vx(ix) > inperi%psizeh(ix)) then
        vx(ix) = vx(ix) - inperi%psize(ix)
        imi(ix) = 0
     else if(vx(ix) < -inperi%psizeh(ix)) then
        vx(ix) = vx(ix) + inperi%psize(ix)
        imi(ix) = 1
     else
        imi(ix) = 2
     end if
  end do

  imirror = 9*imi(1) + 3*imi(2) + imi(3) + 1

end subroutine util_pbneighbor
