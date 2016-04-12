!simu_force_winz
!> @brief Calculates and adds the force related to window potential &
!>        which bind COM of several particles which belong to &
!>        group i to z by a harmonic spring.

subroutine simu_force_winz(irep, force_mp)

  use const_maxsize
  use const_index
  use var_setp, only : inmisc, inwind
  use var_struct, only : xyz_mp_rep, nmp_all, grp
  use var_replica, only : irep2grep, rep2val, inrep
  
  implicit none

  ! ----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: grep, iwind
  integer :: igrp, ilist, imp
  real(PREC) :: z, dx, dy, dz 
  real(PREC) :: kxy, kz
  real(PREC) :: for_x, for_y, for_z
  real(PREC) :: ixyz(3)
  real(PREC) :: coef, dist
  ! ----------------------------------------------------------------------

  grep = irep2grep(irep)
  iwind = inwind%iwinz(grep)

  igrp = inrep%winz_igrp(iwind)
  z    = inrep%winz_z   (iwind)
  kxy  = inrep%winz_kxy (iwind)
  kz   = inrep%winz_kz  (iwind)

  ixyz(:) = 0.0
  do ilist = 1, grp%nmp(igrp)
     ixyz(:) = &
          ixyz(:) + &
          xyz_mp_rep(:, grp%implist(ilist, igrp), irep) * &
          grp%mass_fract(ilist, igrp)
  end do

  dz = ixyz(3) - z

  for_z = -2.0e0_PREC * kz * dz

  do ilist = 1, grp%nmp(igrp)
     imp = grp%implist(ilist, igrp)
     force_mp(3, imp) = force_mp(3, imp) + for_z * grp%mass_fract(ilist, igrp)
  end do

  dx = ixyz(1) - 5500.0e0_PREC
  dy = ixyz(2) - 5500.0e0_PREC

  dist = sqrt(dx * dx + dy * dy)

  coef = -2.0e0_PREC * kxy / dist;

  for_x = coef * dx;
  for_y = coef * dy;

  do ilist = 1, grp%nmp(igrp)
     imp = grp%implist(ilist, igrp)
     force_mp(1, imp) = force_mp(1, imp) + for_x * grp%mass_fract(ilist, igrp)
     force_mp(2, imp) = force_mp(2, imp) + for_y * grp%mass_fract(ilist, igrp)
  end do

end subroutine simu_force_winz
