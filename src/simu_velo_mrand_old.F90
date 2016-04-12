! ***********************************************************************
! Maxwell distribution
! ***********************************************************************
subroutine simu_velo_mrand_old(tempk)

  use const_maxsize, only : PREC, MXMP
  use const_physical, only : BOLTZC
  use var_setp, only : irand, ifix_mp
  use var_struct, only : nmp_real, cmass_mp
  use var_simu, only : velo_mp
  implicit none

  ! --------------------------------------------------------------------
  real(PREC), intent(in) :: tempk

  ! --------------------------------------------------------------------
  ! function
  real(PREC) :: recipe_gasdev

  ! --------------------------------------------------------------------
  ! local variables
  integer :: imp, idimn, irep
  real(PREC) :: coef
  integer, parameter :: MXNOISE = 524288 ! the number of initial call random number


  ! --------------------------------------------------------------------
  ! for the same results
  do imp = 1, MXNOISE
     coef = recipe_gasdev(irand)
  end do

  irep = 1
  do imp = 1, nmp_real
     if(ifix_mp(imp) == 1) then
        velo_mp(1:3, imp, irep) = 0.0
     else
        coef = sqrt(tempk * BOLTZC / cmass_mp(imp))
        do idimn = 1, 3
           velo_mp(idimn, imp, irep) = coef * recipe_gasdev(irand)
        end do
     end if
  end do
!  write (*, *) coef, MXNOISE, nmp_real, irand

end subroutine simu_velo_mrand_old
