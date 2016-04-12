! util_bondangle
!> @brief Calculate the bond angle defined by three particles

subroutine util_bondangle(imp1, imp2, imp3, co_theta, xyz_ba)
      
  use const_maxsize
  implicit none
  
  ! ---------------------------------------------------------------------
  integer, intent(in) :: imp1, imp2, imp3
  real(PREC), intent(out) :: co_theta
  real(PREC), intent(in) :: xyz_ba(3, MXMP)

  ! ---------------------------------------------------------------------
  ! local variables
  real(PREC) :: c11, c21, c22
  real(PREC) :: v21(3), v32(3)
     
  ! ---------------------------------------------------------------------
  v21(1:3) = xyz_ba(1:3, imp2) - xyz_ba(1:3, imp1)
  v32(1:3) = xyz_ba(1:3, imp3) - xyz_ba(1:3, imp2)
  
  c11 = v21(1) * v21(1) + v21(2) * v21(2) + v21(3) * v21(3)  
  c22 = v32(1) * v32(1) + v32(2) * v32(2) + v32(3) * v32(3) 
  c21 = v32(1) * v21(1) + v32(2) * v21(2) + v32(3) * v21(3)
  
  co_theta = -c21 / sqrt(c11 * c22)

  if(co_theta > 1.0e0_PREC) then
     co_theta = 1.0e0_PREC
  else if(co_theta < -1.0e0_PREC) then
     co_theta = -1.0e0_PREC
  end if
  
end subroutine util_bondangle
