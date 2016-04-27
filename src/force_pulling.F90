! force_pulling
!> @brief Calculate force of pulling option

! ************************************************************************
subroutine force_pulling(irep, force_mp)

  use const_maxsize
  use var_setp, only : insimu, inmisc
  use var_struct, only : xyz_mp_rep, nmp_all
  use var_replica, only : irep2grep
  implicit none

  ! ----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: ipull, imp, jmp, grep
!  real(PREC) :: dist
  real(PREC) :: dx, dy, dz, cbd2
  real(PREC) :: force_xyz(1:3)
  
  ! ----------------------------------------------------------------------

  do ipull = 1, inmisc%npull
     force_xyz(1:3) = 0.0e0_PREC

     imp = inmisc%ipull2mp(ipull)

     if(inmisc%coef_pull(ipull) <= 0.0e0_PREC) then
        force_xyz(1:3) = inmisc%pull_xyz(1:3, ipull)
     else
        inmisc%pu_xyz(1:3, ipull) = inmisc%pu_xyz(1:3, ipull) + &
             inmisc%pull_xyz(1:3, ipull) * insimu%tstep_size

        dx = xyz_mp_rep(1, imp, irep) - inmisc%pu_xyz(1, ipull)
        dy = xyz_mp_rep(2, imp, irep) - inmisc%pu_xyz(2, ipull)
        dz = xyz_mp_rep(3, imp, irep) - inmisc%pu_xyz(3, ipull)
!        dist = (dx**2 + dy**2 + dz**2)**0.5

        cbd2 = -2.0e0_PREC * inmisc%coef_pull(ipull)

        force_xyz(1) = cbd2 * dx
        force_xyz(2) = cbd2 * dy
        force_xyz(3) = cbd2 * dz
     end if

     force_mp(1:3, imp) = force_mp(1:3, imp) + force_xyz(1:3)
  end do


  do ipull = 1, inmisc%npull_unravel
     grep = irep2grep(irep)
     force_xyz(:) = inmisc%pull_unravel_xyz(:,ipull,grep)
     imp = inmisc%ipull_unravel2mp(1, ipull)
     jmp = inmisc%ipull_unravel2mp(2, ipull)
     force_mp(:, imp) = force_mp(:, imp) + force_xyz(:)
     force_mp(:, jmp) = force_mp(:, jmp) - force_xyz(:)
  enddo

end subroutine force_pulling
