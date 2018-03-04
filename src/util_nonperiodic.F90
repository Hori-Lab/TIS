! util_periodic
!> @brief back mapping from pxyz_mp_rep to xyz_mp_rep which is used for local
! interations such as bond

subroutine util_nonperiodic(irep)
      
  use const_maxsize
  use var_setp,   only : inperi
  use var_struct, only : nmp_real, pxyz_mp_rep, xyz_mp_rep
  implicit none
  
  ! ---------------------------------------------------------------------
  integer, intent(in) :: irep

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: imp, ix
  real(PREC) :: pwide(3), pmax(3), pmin(3)
  real(PREC) :: dx(3)
  real(PREC) :: tx, d
     
  ! ---------------------------------------------------------------------
  
  pwide(1:3) = inperi%psize(1:3)
  pmax(1:3) = 0.5*pwide(1:3)
  pmin(1:3) = -0.5*pwide(1:3)

  dx(1:3) = 0.0
  xyz_mp_rep(:,1,irep) = pxyz_mp_rep(:,1,irep)
  
  do imp = 2, nmp_real
     do ix = 1, 3
        tx = pxyz_mp_rep(ix, imp, irep) + dx(ix)
        d = tx - xyz_mp_rep(ix,imp-1,irep)
        do while (.true.)
           if(d > pmax(ix)) then
              dx(ix) = dx(ix) - inperi%psize(ix)
           else if(d < pmin(ix)) then
              dx(ix) = dx(ix) + inperi%psize(ix)
           else
              exit
           end if
           tx = pxyz_mp_rep(ix, imp, irep) + dx(ix)
           d = tx - xyz_mp_rep(ix,imp-1,irep)
        enddo
        xyz_mp_rep(ix,imp,irep) = tx
     end do
  end do

end subroutine util_nonperiodic
