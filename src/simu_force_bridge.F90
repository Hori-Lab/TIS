!simu_force_bridge
!> @brief Calculates and adds the force related to bridge interactions &
!>        which bind two mass-points by a harmonic spring.

subroutine simu_force_bridge(irep, force_mp)

  use const_maxsize
  use var_setp, only : inmisc
  use var_struct, only : xyz_mp_rep, nmp_all, grp, cmass_mp
  use var_simu, only : flg_ppr_release
  implicit none

  ! ----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: i, ibrid, imp, jmp, igrp, jgrp
  real(PREC) :: d(3), dist, cbd2, for
  real(PREC) :: force(3)
  real(PREC) :: ixyz(3), jxyz(3)
  
  ! ----------------------------------------------------------------------
  ! BRIDGE
  do ibrid = 1, inmisc%nbrid
     force(:) = 0.0e0_PREC

     imp = inmisc%ibrid2mp(1, ibrid)
     jmp = inmisc%ibrid2mp(2, ibrid)
     d(:) = xyz_mp_rep(:, imp, irep) - xyz_mp_rep(:, jmp, irep)

     dist = sqrt(dot_product(d,d))
     if((inmisc%i_lower_bound == 0 .and. dist > inmisc%brid_dist(ibrid)) &
          .or. inmisc%i_lower_bound /= 0) then
        cbd2 = -2.0e0_PREC * inmisc%coef_brid(ibrid)
        for = cbd2 * (dist - inmisc%brid_dist(ibrid)) / dist
        force(:) = for * d(:)
     else
        exit
     end if

     force_mp(:, imp) = force_mp(:, imp) + force(:)
     force_mp(:, jmp) = force_mp(:, jmp) - force(:)
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

     if (inmisc%i_lower_bound /= 0 .or. &
         (inmisc%i_lower_bound == 0 .and. dist > inmisc%brid_com_dist(ibrid))) then
        cbd2 = -2.0e0_PREC * inmisc%coef_brid_com(ibrid)
        for = cbd2 * (dist - inmisc%brid_com_dist(ibrid)) / dist
        force(:) = for * d(:)
     else
        exit
     end if

     do i = 1, grp%nmp(igrp)
        imp = grp%implist(i,igrp)
        force_mp(:, imp) = force_mp(:, imp) + force(:) * grp%mass_fract(i,igrp)
     enddo
     do i = 1, grp%nmp(jgrp)
        jmp = grp%implist(i,jgrp)
        force_mp(:, jmp) = force_mp(:, jmp) - force(:) * grp%mass_fract(i,jgrp)
     enddo
  end do

end subroutine simu_force_bridge

