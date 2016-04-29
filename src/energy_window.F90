! energy_window
!> @brief Calculate energy of window option

! ************************************************************************
subroutine energy_window(irep, energy_unit, energy)

  use const_maxsize
  use const_index
  use var_setp,    only : inwind
  use var_struct,  only : xyz_mp_rep, imp2unit, nunit_all
  use var_replica, only : irep2grep, rep2val, inrep
  

  implicit none
  ! ------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy_unit(nunit_all, nunit_all, E_TYPE%MAX)
  real(PREC), intent(inout) :: energy(E_TYPE%MAX)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: grep, iwind
  integer :: imp, jmp, iunit, junit
  real(PREC) :: dx, dy, dz, dist, efull
  real(PREC) :: coef, length
  
  ! ----------------------------------------------------------------------

  ! Obtain window id
  grep = irep2grep(irep)
  iwind = inwind%iwind(grep)

  ! Refer parameters from window id
  coef = inrep%window_property(iwind, WINDTYPE%COEF)
  length = inrep%window_property(iwind, WINDTYPE%LENGTH)
  imp = inrep%window_mp_id(iwind, WINDTYPE%IMP)
  jmp = inrep%window_mp_id(iwind, WINDTYPE%JMP)

  ! Calculate distance between imp and jmp
  dx = xyz_mp_rep(1, imp, irep) - xyz_mp_rep(1, jmp, irep)
  dy = xyz_mp_rep(2, imp, irep) - xyz_mp_rep(2, jmp, irep)
  dz = xyz_mp_rep(3, imp, irep) - xyz_mp_rep(3, jmp, irep)
   
  dist = sqrt(dx**2 + dy**2 + dz**2)

  ! Calculate energy
  efull = coef * (dist - length)**2

  ! Increment energy
  energy(E_TYPE%WINDOW) = energy(E_TYPE%WINDOW) + efull

  iunit = imp2unit(imp)
  junit = imp2unit(jmp)

  energy_unit(iunit, junit, E_TYPE%WINDOW) =   &
       energy_unit(iunit, junit, E_TYPE%WINDOW) + efull
  
end subroutine energy_window
