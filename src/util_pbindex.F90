! util_pbindex
!> @brief calc the neighbor index in periodic boundary condition

subroutine util_pbindex(ineigh2mp, nneigh, irep)

  use const_maxsize
  use var_setp,   only : inperi
  use var_struct, only : pxyz_mp_rep

  implicit none

  ! --------------------------------------------------------------------
  integer, intent(inout) :: ineigh2mp(3, nneigh)
  integer, intent(inout) :: nneigh
  integer, intent(in) :: irep

  ! --------------------------------------------------------------------
  ! local variables
  integer :: imp1, imp2, ineigh
  integer :: ix, imirror
  integer :: imi(3)
  real(PREC) :: vx(3)

  ! --------------------------------------------------------------------

  do ineigh = 1, nneigh
     imp1 = ineigh2mp(1, ineigh)
     imp2 = ineigh2mp(2, ineigh)
 
     vx(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)

     do ix = 1, 3
        if(vx(ix) > inperi%psizeh(ix)) then
!           vx(ix) = vx(ix) - inperi%psize(ix)
           imi(ix) = 1
        else if(vx(ix) < -inperi%psizeh(ix)) then
!           vx(ix) = vx(ix) + inperi%psize(ix)
           imi(ix) = 2
        else
           imi(ix) = 0
        end if
     end do
     imirror = 9*imi(1) + 3*imi(2) + imi(3)

     ineigh2mp(3, ineigh) = imirror
  end do


end subroutine util_pbindex
