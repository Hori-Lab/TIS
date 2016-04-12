! ***********************************************************************
! create a random number distributed uniformly in (0,1).
! The 'ran2' routine in Numerical Recipes-2nd edition.
! ***********************************************************************
! function ran2(idum)
function recipe_rand01(idum)

  use const_maxsize, only : PREC
  implicit none

  real(PREC) :: recipe_rand01

  ! --------------------------------------------------------------------
  integer, intent(inout) :: idum

  ! --------------------------------------------------------------------
  ! local variables
  integer, parameter :: IM1 = 2147483563
  integer, parameter :: IM2 = 2147483399
  integer, parameter :: IMM1 = IM1 - 1
  integer, parameter :: IA1 = 40014
  integer, parameter :: IA2 = 40692
  integer, parameter :: IQ1 = 53668
  integer, parameter :: IQ2 = 52774
  integer, parameter :: IR1 = 12211
  integer, parameter :: IR2 = 3791
  integer, parameter :: NTAB = 32
  integer, parameter :: NDIV = 1 + IMM1 / NTAB
  real(PREC), parameter :: AM = 1.0e0_PREC / IM1
  real(PREC), parameter :: EPS = 1.2e-7_PREC
  real(PREC), parameter :: RNMX = 1.0e0_PREC - EPS


  integer :: j, k
  integer, save :: iy = 0
  integer, save :: idum2 = 123456789
  integer, save :: iv(NTAB)
  data iv / NTAB*0 /
  real(PREC) :: ran2

  ! --------------------------------------------------------------------
  if (idum <= 0) then
     idum = max(-idum, 1)
     idum2 = idum
     do j = NTAB + 8, 1, -1
        k = idum / IQ1
        idum = IA1 * (idum - k * IQ1) - k * IR1
        if (idum < 0) idum = idum + IM1
        if (j <= NTAB) iv(j) = idum
     end do
     iy = iv(1)
  end if

  k = idum / IQ1
  idum = IA1 * (idum - k * IQ1) - k * IR1

  if (idum < 0) idum = idum + IM1
  k = idum2 / IQ2
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2

  if (idum2 < 0) idum2 = idum2 + IM2
  j = 1 + iy / NDIV
  iy = iv(j) - idum2
  iv(j) = idum

  if(iy < 1) iy = iy + IMM1
  ran2 = min(AM * iy, RNMX)
  recipe_rand01 = ran2

end function recipe_rand01
  
