! energy_pulling
!> @brief Calculate energy of pulling option

! ************************************************************************
subroutine energy_pulling(irep, energy_unit, energy)

  use const_maxsize
  use const_index
  use var_setp,    only : inmisc
  use var_struct,  only : xyz_mp_rep, imp2unit
  use var_replica, only : irep2grep

  implicit none
  ! ----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(inout) :: energy_unit(:,:,:) ! (unit, unit, E_TYPE%MAX)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: ipull, imp, jmp, iunit, junit, grep
  real(PREC) :: v(3), cbd2, efull, force_xyz(3)
  
  ! ----------------------------------------------------------------------
  do ipull = 1, inmisc%npull
   
     imp = inmisc%ipull2mp(ipull)

     if(inmisc%coef_pull(ipull) <= 0.0e0_PREC) then
   
     else
        v(:) = xyz_mp_rep(:, imp, irep) - inmisc%pu_xyz(:, ipull)

        efull = inmisc%coef_pull(ipull) * dot_product(v,v)
   
        energy(E_TYPE%PULLING) = energy(E_TYPE%PULLING) + efull
        iunit = imp2unit(imp)
        junit = iunit
        energy_unit(iunit, junit, E_TYPE%PULLING) = energy_unit(iunit, junit, E_TYPE%PULLING) + efull
     end if

  enddo

  do ipull = 1, inmisc%npull_unravel
     imp = inmisc%ipull_unravel2mp(1, ipull)
     jmp = inmisc%ipull_unravel2mp(2, ipull)

     v(:) = xyz_mp_rep(:,imp,irep) - xyz_mp_rep(:,jmp,irep)
     grep = irep2grep(irep)
     force_xyz(:) = inmisc%pull_unravel_xyz(:,ipull,grep)
     efull = - dot_product(v, force_xyz)

     energy(E_TYPE%PULLING) = energy(E_TYPE%PULLING) + efull
     iunit = imp2unit(imp)
     junit = imp2unit(jmp)
     energy_unit(iunit, junit, E_TYPE%PULLING) = energy_unit(iunit, junit, E_TYPE%PULLING) + 0.5 * efull
  enddo

end subroutine energy_pulling
