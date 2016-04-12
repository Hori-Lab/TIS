!rfunc_cubiceq
!> @brief  Solves a cubic equation by Cardano-Tartaglia method  &
!>         and return the smallest real root.

!***********************************************************************
!                                                                      
!   Solve a cubic equation by Cardano-Tartaglia method
!            and return the smallest real root:
!          
!             a*x**3 + b*x**2 + c*x + d = 0
!
!***********************************************************************

function rfunc_cubiceq(a, b, c, d)

  use const_maxsize
  use const_physical
  implicit none

  ! -----------------------------------------------------
  real(PREC) :: rfunc_cubiceq
  real(PREC), intent(in) :: a, b, c, d

  !------------------------------------------------------
  real(PREC) :: ca, cb, cc, p, q, disc, &
                s1, s2, r1, r2, r3, rmin, &
                phi, s11, s22

#define IND  3.0e0_PREC
#define IND2 9.0e0_PREC
#define IND3 27.0e0_PREC
!------------------------------------------------------
  
  ca = b / a
  cb = c / a
  cc = d / a
  
  p = cb - (ca**2) / IND
  q = cc + (2.0e0_PREC * ca**3 - IND2 * ca * cb) / IND3
  
  ! Discriminant
  disc = 0.25e0_PREC * q**2 + (p**3) / IND3
  if(disc >= 0.0e0_PREC) then
    s11 = -0.5e0_PREC * q + sqrt(disc)
    s22 = -0.5e0_PREC * q - sqrt(disc)
    s1 = abs(s11)**(1.0e0_PREC / IND)
    s2 = abs(s22)**(1.0e0_PREC / IND)
    if(s11 < 0) s1 = -s1
    if(s22 < 0) s2 = -s2
   
    ! D>0, one real root and two conjugated complex roots
    rmin = s1 + s2

    ! D=0, three real roots and at least two are equal
    if(disc <= ZERO_JUDGE) then
      r2 = -s1
      if(r2 < rmin) rmin = r2
    end if

  else
    ! D<0, three real roots
    s1 = sqrt(abs(p)**3 / IND3)
    s2 = sqrt(abs(p) / IND) 
    phi   = acos(-0.5e0_PREC * q / s1)
    
    r1 =  2.0e0_PREC * s2 * cos(phi / IND)
    r2 = -2.0e0_PREC * s2 * cos((phi + F_PI) / IND)
    r3 = -2.0e0_PREC * s2 * cos((phi - F_PI) / IND)
    
    rmin = r1
    if(r2 < rmin) rmin = r2
    if(r3 < rmin) rmin = r3
  end if

  rfunc_cubiceq = rmin - ca / IND

end function rfunc_cubiceq
