! simu_energy_bridge
!> @brief Calculate energy of bridge option

! ************************************************************************
subroutine simu_energy_bridge(irep, pnle_unit, pnlet)

  use const_maxsize
  use const_index
  use var_setp,    only : inmisc
  use var_struct,  only : xyz_mp_rep, imp2unit, grp
  use var_simu,    only : flg_ppr_release

  implicit none
  ! ------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: pnle_unit(:,:,:) ! (unit, unit, E_TYPE%MAX)
  real(PREC), intent(inout) :: pnlet(:)         ! (E_TYPE%MAX)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: i, ibrid, imp, jmp, iunit, junit
  integer :: igrp, jgrp
  real(PREC) :: dx, dy, dz, dist, cbd2, efull
  real(PREC) :: ixyz(3), jxyz(3), d(3)
  
  ! ----------------------------------------------------------------------
  ! BRIDGE
  do ibrid = 1, inmisc%nbrid

     imp = inmisc%ibrid2mp(1, ibrid)
     jmp = inmisc%ibrid2mp(2, ibrid)
     dx = xyz_mp_rep(1, imp, irep) - xyz_mp_rep(1, jmp, irep)
     dy = xyz_mp_rep(2, imp, irep) - xyz_mp_rep(2, jmp, irep)
     dz = xyz_mp_rep(3, imp, irep) - xyz_mp_rep(3, jmp, irep)
   
     dist = sqrt(dx**2 + dy**2 + dz**2)

     if(     (inmisc%i_lower_bound == 0 .and. dist > inmisc%brid_dist(ibrid)) &
         .or. inmisc%i_lower_bound == 1 &
         .or.(inmisc%i_lower_bound == 2 .and. dist < inmisc%brid_dist(ibrid)) ) then
        cbd2 = inmisc%coef_brid(ibrid)
        efull = cbd2 * (dist - inmisc%brid_dist(ibrid))**2
        
        pnlet(E_TYPE%BRIDGE) = pnlet(E_TYPE%BRIDGE) + efull
        iunit = imp2unit(imp)
        junit = imp2unit(jmp)
        pnle_unit(iunit, junit, E_TYPE%BRIDGE) =   &
              pnle_unit(iunit, junit, E_TYPE%BRIDGE) + efull
     end if
   
  end do
     
  ! BRIDGE_CENTER
  do ibrid = 1, inmisc%nbrid_com

     if (flg_ppr_release(ibrid)) then
        cycle
     endif

     igrp = inmisc%ibrid_com2grp(1, ibrid)
     jgrp = inmisc%ibrid_com2grp(2, ibrid)

     ixyz(:) = 0.0
     do i = 1, grp%nmp(igrp)
        ixyz(:) = ixyz(:) + xyz_mp_rep(:, grp%implist(i,igrp), irep) * grp%mass_fract(i,igrp)
     enddo

     jxyz(:) = 0.0
     do i = 1, grp%nmp(jgrp)
        jxyz(:) = jxyz(:) + xyz_mp_rep(:, grp%implist(i,jgrp), irep) * grp%mass_fract(i,jgrp)
     enddo

     d(:) = ixyz(:) - jxyz(:)
     dist = sqrt(dot_product(d,d))

     if((inmisc%i_lower_bound == 0 .and. dist > inmisc%brid_com_dist(ibrid)) &
          .or. inmisc%i_lower_bound /= 0) then
        cbd2 = inmisc%coef_brid_com(ibrid)
        efull = cbd2 * (dist - inmisc%brid_com_dist(ibrid))**2
        
        pnlet(E_TYPE%BRIDGE) = pnlet(E_TYPE%BRIDGE) + efull

        !! NOTE: Energy for unit is not considered for BRIDGE_CENTER.

        !iunit = imp2unit(imp)
        !junit = imp2unit(jmp)
        !pnle_unit(iunit, junit, E_TYPE%BRIDGE) =   &
        !      pnle_unit(iunit, junit, E_TYPE%BRIDGE) + efull
     end if
   
  end do
end subroutine simu_energy_bridge
