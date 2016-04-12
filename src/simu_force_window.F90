!simu_force_window
!> @brief Calculates and adds the force related to window potential &
!>        which bind two mass-points by a harmonic spring.

subroutine simu_force_window(irep, force_mp)

  use const_maxsize
  use const_index
  use var_setp, only : inmisc, inwind
  use var_struct, only : xyz_mp_rep, nmp_all
  use var_replica, only : irep2grep, rep2val, inrep
  implicit none

  ! ----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: grep
  integer :: imp, jmp
  integer :: iwind
  real(PREC) :: dx, dy, dz, dist, cbd2, for
  real(PREC) :: force_x, force_y, force_z
  real(PREC) :: coef, length
  
  ! ----------------------------------------------------------------------

  grep = irep2grep(irep)
  iwind = inwind%iwind(grep)

  coef = inrep%window_property(iwind, WINDTYPE%COEF)
  length = inrep%window_property(iwind, WINDTYPE%LENGTH)
  imp = inrep%window_mp_id(iwind, WINDTYPE%IMP)
  jmp = inrep%window_mp_id(iwind, WINDTYPE%JMP)

  dx = xyz_mp_rep(1, imp, irep) - xyz_mp_rep(1, jmp, irep)
  dy = xyz_mp_rep(2, imp, irep) - xyz_mp_rep(2, jmp, irep)
  dz = xyz_mp_rep(3, imp, irep) - xyz_mp_rep(3, jmp, irep)
  
  dist = sqrt(dx**2 + dy**2 + dz**2)
  
  cbd2 = -2.0e0_PREC * coef
  for = cbd2 * (dist - length) / dist
  
  force_x = for * dx
  force_y = for * dy
  force_z = for * dz 
    
  force_mp(1, imp) = force_mp(1, imp) + force_x
  force_mp(2, imp) = force_mp(2, imp) + force_y 
  force_mp(3, imp) = force_mp(3, imp) + force_z
  force_mp(1, jmp) = force_mp(1, jmp) - force_x
  force_mp(2, jmp) = force_mp(2, jmp) - force_y 
  force_mp(3, jmp) = force_mp(3, jmp) - force_z

end subroutine simu_force_window

