!force_window
!> @brief Calculates and adds the force related to window potential &
!>        which bind two mass-points by a harmonic spring.

subroutine force_window(irep, force_mp)

  use const_maxsize
  use const_index
  use var_setp, only : inwind
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
  real(PREC) :: v(3), dist, for(3)
  real(PREC) :: coef, length
  
  ! ----------------------------------------------------------------------

  grep = irep2grep(irep)
  iwind = inwind%iwind(grep)

  coef = inrep%window_property(iwind, WINDTYPE%COEF)
  length = inrep%window_property(iwind, WINDTYPE%LENGTH)
  imp = inrep%window_mp_id(iwind, WINDTYPE%IMP)
  jmp = inrep%window_mp_id(iwind, WINDTYPE%JMP)

  v(:) = xyz_mp_rep(:, imp, irep) - xyz_mp_rep(:, jmp, irep)
  
  dist = norm2(v)
  
  for(:) = -2.0e0_PREC * coef * (dist - length) / dist * v(:)
    
  force_mp(:, imp) = force_mp(:, imp) + for(:)
  force_mp(:, jmp) = force_mp(:, jmp) - for(:)

end subroutine force_window

