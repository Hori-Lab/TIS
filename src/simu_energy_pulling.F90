! simu_energy_pulling
!> @brief Calculate energy of pulling option

! ************************************************************************
subroutine simu_energy_pulling(irep, e_exv_unit, e_exv)

  use const_maxsize
  use const_index
  use var_setp,    only : inmisc
  use var_struct,  only : xyz_mp_rep, imp2unit
  use var_replica, only : irep2grep

  implicit none
  ! ----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: e_exv(:)         ! (E_TYPE%MAX)
  real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (unit, unit, E_TYPE%MAX)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: ipull, imp, jmp, iunit, junit, grep
  real(PREC) :: dx, dy, dz, dist, cbd2, efull, vij(3), force_xyz(3)
  
  ! ----------------------------------------------------------------------
  do ipull = 1, inmisc%npull
   
     imp = inmisc%ipull2mp(ipull)

     if(inmisc%coef_pull(ipull) <= 0.0e0_PREC) then
   
     else
        dx = xyz_mp_rep(1, imp, irep) - inmisc%pu_xyz(1, ipull)
        dy = xyz_mp_rep(2, imp, irep) - inmisc%pu_xyz(2, ipull)
        dz = xyz_mp_rep(3, imp, irep) - inmisc%pu_xyz(3, ipull)
        dist = sqrt(dx**2 + dy**2 + dz**2)
   
        cbd2 = inmisc%coef_pull(ipull)
        efull = cbd2 * dist**2;
   
        e_exv(E_TYPE%PULLING) = e_exv(E_TYPE%PULLING) + efull
        iunit = imp2unit(imp)
        junit = iunit
        e_exv_unit(iunit, junit, E_TYPE%PULLING) = e_exv_unit(iunit, junit, E_TYPE%PULLING) + efull
     end if

  enddo

  do ipull = 1, inmisc%npull_unravel
     imp = inmisc%ipull_unravel2mp(1, ipull)
     jmp = inmisc%ipull_unravel2mp(2, ipull)

     vij(:) = xyz_mp_rep(:,imp,irep) - xyz_mp_rep(:,jmp,irep)
     grep = irep2grep(irep)
     force_xyz(:) = inmisc%pull_unravel_xyz(:,ipull,grep)
     efull = - dot_product(vij, force_xyz)

     e_exv(E_TYPE%PULLING) = e_exv(E_TYPE%PULLING) + efull
     iunit = imp2unit(imp)
     junit = imp2unit(jmp)
     e_exv_unit(iunit, junit, E_TYPE%PULLING) = e_exv_unit(iunit, junit, E_TYPE%PULLING) + 0.5 * efull
  enddo

end subroutine simu_energy_pulling
