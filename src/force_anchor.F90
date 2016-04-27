!force_anchor
!> @brief Calculates and adds the force related to constrain by anchor.

subroutine force_anchor(irep, force_mp)

  use const_maxsize
  use var_setp, only : inmisc
  use var_struct, only : xyz_mp_rep, nmp_all, grp
  implicit none

  ! ----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: ianc, ilist, imp, igrp
  real(PREC) :: dx, dy, dz, dist, cbd2, for
  real(PREC) :: force_x, force_y, force_z
  real(PREC) :: ixyz(3)
  ! ----------------------------------------------------------------------
  
  do ianc = 1, inmisc%nanc
     force_x = 0.0e0_PREC
     force_y = 0.0e0_PREC
     force_z = 0.0e0_PREC

     imp = inmisc%ianc2mp(ianc)
     dx = xyz_mp_rep(1, imp, irep) - inmisc%anc_xyz(1, ianc)
     dy = xyz_mp_rep(2, imp, irep) - inmisc%anc_xyz(2, ianc)
     dz = xyz_mp_rep(3, imp, irep) - inmisc%anc_xyz(3, ianc)

     dist = sqrt(dx**2 + dy**2 + dz**2)
     if(dist > inmisc%anc_dist(ianc)) then
        cbd2 = -2.0e0_PREC * inmisc%coef_anc(ianc)
        for = cbd2 * (dist - inmisc%anc_dist(ianc)) / dist
        force_x = for * dx
        force_y = for * dy
        force_z = for * dz
     end if

     force_mp(1, imp) = force_mp(1, imp) + force_x
     force_mp(2, imp) = force_mp(2, imp) + force_y 
     force_mp(3, imp) = force_mp(3, imp) + force_z
  end do

  do ianc = 1, inmisc%nanc_com_ini
     
     force_x = 0.0e0_PREC
     force_y = 0.0e0_PREC
     force_z = 0.0e0_PREC
     
     igrp = inmisc%ianc_com_ini2grp(ianc)
     
     ixyz(:) = 0.0
     do ilist = 1, grp%nmp(igrp)
        ixyz(:) = &
             ixyz(:) + &
             xyz_mp_rep(:, grp%implist(ilist, igrp), irep) * &
             grp%mass_fract(ilist, igrp)
     end do

     dx = ixyz(1) - inmisc%anc_com_ini_xyz(1, ianc, irep)
     dy = ixyz(2) - inmisc%anc_com_ini_xyz(2, ianc, irep)
     dz = ixyz(3) - inmisc%anc_com_ini_xyz(3, ianc, irep)

     dist = sqrt(dx**2 + dy**2 + dz**2)

     if(dist > inmisc%anc_com_ini_dist(ianc)) then
        cbd2 = -2.0e0_PREC * inmisc%coef_anc_com_ini(ianc)
        for = cbd2 * (dist - inmisc%anc_com_ini_dist(ianc)) / dist
        force_x = for * dx
        force_y = for * dy
        force_z = for * dz
     end if

     do ilist = 1, grp%nmp(igrp)
        imp = grp%implist(ilist, igrp)
        force_mp(1, imp) = &
             force_mp(1, imp) + force_x * grp%mass_fract(ilist, igrp)
        force_mp(2, imp) = &
             force_mp(2, imp) + force_y * grp%mass_fract(ilist, igrp)
        force_mp(3, imp) = &
             force_mp(3, imp) + force_z * grp%mass_fract(ilist, igrp)
     end do

  end do

end subroutine force_anchor
