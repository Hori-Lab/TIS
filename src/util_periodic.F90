! util_periodic
!> @brief apply periodic boundary

subroutine util_periodic(irep)
      
  use const_maxsize
  use var_inp, only : inperi
  use var_struct, only : nmp_real, pxyz_mp_rep
  implicit none
  
  ! ---------------------------------------------------------------------
  integer, intent(in) :: irep

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: imp, ix
  real(PREC) :: pwide(3), pmax(3), pmin(3)
  real(PREC) :: tx
     
  ! ---------------------------------------------------------------------
  
  pwide(1:3) = inperi%psize(1:3)
  pmax(1:3) = 0.5*pwide(1:3)
  pmin(1:3) = -0.5*pwide(1:3)
  
  do imp = 1, nmp_real
     do ix = 1, 3
        tx = pxyz_mp_rep(ix, imp, irep)
        if(tx > pmax(ix)) then
           tx = tx - pwide(ix) * (int((tx - pmax(ix))/pwide(ix)) + 1)
        else if(tx < pmin(ix)) then
           tx = tx + pwide(ix) * (int((pmin(ix) - tx)/pwide(ix)) + 1)
        end if
        pxyz_mp_rep(ix, imp, irep) = tx
     end do
  end do
  

end subroutine util_periodic
