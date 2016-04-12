! ***********************************************************************
! function: real(PREC) recipe_gasdev(integer)
! make Gaussian white noise  (average=0, variance=1)
! ***********************************************************************
function recipe_gasdev(iseed)

  use const_maxsize, only : PREC
  implicit none

  ! --------------------------------------------------------------------
  real(PREC) :: recipe_gasdev

  ! --------------------------------------------------------------------
  integer, intent(inout) :: iseed

  ! --------------------------------------------------------------------
  ! function
  real(PREC) :: recipe_rand01

  ! --------------------------------------------------------------------
  ! local variables
  integer, save :: iset = 0
  real(PREC) :: fac, r, v1, v2
  real(PREC), save :: gset = 0.0e0_PREC
  
  ! --------------------------------------------------------------------
  recipe_gasdev = 0.0e0_PREC

  if(iset == 0) then
     do
        v1 = 2.0e0_PREC * recipe_rand01(iseed) - 1.0e0_PREC
        v2 = 2.0e0_PREC * recipe_rand01(iseed) - 1.0e0_PREC
        r = v1**2 + v2**2
        if(r < 1.0e0_PREC .and. r > 0.0e0_PREC) exit
     end do
     fac = sqrt(-2.0e0_PREC * log(r) / r)
     gset = v1 * fac
     iset = 1
     recipe_gasdev = v2 * fac

  else
     iset = 0
     recipe_gasdev = gset
  end if
  
end function recipe_gasdev
