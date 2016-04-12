! util_dihangle
!> @brief This subroutine is to calculate the dihedral angle from the coordinate data.

! ************************************************************************
subroutine util_dihangle(imp1, imp2, imp3, imp4, dih_angle, co_dih, si_dih, xyz_dih)
  
  use const_maxsize
  use const_physical
  implicit none
  
  ! ---------------------------------------------------------------------
  integer, intent(in) :: imp1, imp2, imp3, imp4
  real(PREC), intent(out) :: dih_angle, co_dih, si_dih
  real(PREC), intent(in) :: xyz_dih(3, MXMP)

  ! ---------------------------------------------------------------------
  ! local variables
!  real(PREC) :: pi2
  real(PREC) :: c11, c12, c13, c22, c23, c33
  real(PREC) :: t1, t3, t4, t3t4, zahyokei
  real(PREC) :: v21(3), v32(3), v43(3)

  ! ---------------------------------------------------------------------
!  pi2 = 2.0e0_PREC * F_PI

  v21(1:3) = xyz_dih(1:3, imp2) - xyz_dih(1:3, imp1)
  v32(1:3) = xyz_dih(1:3, imp3) - xyz_dih(1:3, imp2)
  v43(1:3) = xyz_dih(1:3, imp4) - xyz_dih(1:3, imp3)

  c11 = v21(1) * v21(1) + v21(2) * v21(2) + v21(3) * v21(3)
  c22 = v32(1) * v32(1) + v32(2) * v32(2) + v32(3) * v32(3)
  c33 = v43(1) * v43(1) + v43(2) * v43(2) + v43(3) * v43(3)
  
  c12 = v21(1) * v32(1) + v21(2) * v32(2) + v21(3) * v32(3)
  c13 = v21(1) * v43(1) + v21(2) * v43(2) + v21(3) * v43(3)
  c23 = v32(1) * v43(1) + v32(2) * v43(2) + v32(3) * v43(3)
  
  t1 = c12 * c23 - c13 * c22
  t3 = c11 * c22 - c12 * c12
  t4 = c22 * c33 - c23 * c23

  if (t3 < ZERO_JUDGE) t3 = ZERO_JUDGE
  if (t4 < ZERO_JUDGE) t4 = ZERO_JUDGE

  t3t4 = sqrt(t3 * t4)
  co_dih = t1 / t3t4
  
  if(co_dih > 1.0e0_PREC) then ! when co > 1 ,or < 1
     co_dih = 1.0e0_PREC
  else if(co_dih < -1.0e0_PREC) then
     co_dih = -1.0e0_PREC
  end if
  
  zahyokei = v21(1) * v32(2) * v43(3) + &
       v32(1) * v43(2) * v21(3) + &
       v43(1) * v21(2) * v32(3) - &
       v43(1) * v32(2) * v21(3) - &
       v21(1) * v43(2) * v32(3) - &
       v32(1) * v21(2) * v43(3)

  si_dih = sqrt(1.0e0_PREC - co_dih**2)
  if (si_dih < ZERO_JUDGE) si_dih = ZERO_JUDGE

  if(zahyokei < 0.0e0_PREC) then
     si_dih = - si_dih
  end if

  if(zahyokei > 0.0e0_PREC) then
     dih_angle = acos(co_dih)
  else
     dih_angle = - acos(co_dih)
  end if
  
end subroutine util_dihangle
