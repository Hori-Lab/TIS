! simu_energy_anchor
!> @brief Calculate energy of anchor option

! ************************************************************************
subroutine simu_energy_anchor(irep, pnle_unit, pnlet)

  use const_maxsize
  use const_index
  use var_setp, only : inmisc
  use var_struct, only : xyz_mp_rep, imp2unit, grp
  implicit none

  ! ------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: pnle_unit(:,:,:) ! (unit, unit, E_TYPE%MAX)
  real(PREC), intent(inout) :: pnlet(:)         ! (E_TYPE%MAX)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: ianc, imp, iunit, junit, igrp, ilist
  real(PREC) :: dx, dy, dz, dist, cbd2, efull
  real(PREC) :: ixyz(3)
  ! ----------------------------------------------------------------------
  do ianc = 1, inmisc%nanc
   
     imp = inmisc%ianc2mp(ianc)
     dx = xyz_mp_rep(1, imp, irep) - inmisc%anc_xyz(1, ianc)
     dy = xyz_mp_rep(2, imp, irep) - inmisc%anc_xyz(2, ianc)
     dz = xyz_mp_rep(3, imp, irep) - inmisc%anc_xyz(3, ianc)

     dist = sqrt(dx**2 + dy**2 + dz**2)
     if(dist > inmisc%anc_dist(ianc)) then
        cbd2 = inmisc%coef_anc(ianc)
        efull = cbd2 * (dist - inmisc%anc_dist(ianc))**2
   
        pnlet(E_TYPE%ANCHOR) = pnlet(E_TYPE%ANCHOR) + efull
        iunit = imp2unit(imp)
        junit = iunit
        pnle_unit(iunit, junit, E_TYPE%ANCHOR) =  &
                   pnle_unit(iunit, junit, E_TYPE%ANCHOR) + efull
     end if

  enddo

  do ianc = 1, inmisc%nanc_com_ini

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
        cbd2 = inmisc%coef_anc_com_ini(ianc)
        efull = cbd2 * (dist - inmisc%anc_com_ini_dist(ianc))**2
        pnlet(E_TYPE%ANCHOR) = pnlet(E_TYPE%ANCHOR) + efull
     end if

  end do
  
end subroutine simu_energy_anchor
