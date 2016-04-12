! simu_velo_settemp
!> @brief Control temperature for Brendsen thermostat

! ***********************************************************************
subroutine simu_velo_settemp(velo_mp, irep, tempk)

  use const_maxsize
  use const_physical
  use var_setp, only : inpara, insimu, ifix_mp
  use var_struct, only : nmp_real, cmass_mp
  implicit none
      
  ! --------------------------------------------------------------------
  real(PREC), intent(in) :: tempk
  real(PREC), intent(inout) :: velo_mp(:,:,:)
  integer, intent(in) :: irep

  ! --------------------------------------------------------------------
  ! local variables
  integer :: imp, kmp
  integer :: njiyu
  real(PREC) :: ttemp, temp_now, veloet, velo_adjst
  real(PREC) :: rtarget

  ! --------------------------------------------------------------------
  ttemp = tempk * BOLTZC
  velo_adjst = inpara%velo_adjst

  veloet = 0.0
  kmp = 0
  do imp = 1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     veloet = veloet + 0.5e0_PREC * cmass_mp(imp) * (velo_mp(1, imp, irep)**2 &
          + velo_mp(2, imp, irep)**2 + velo_mp(3, imp, irep)**2)
     kmp = kmp + 1
  end do
  
  if(insimu%i_no_trans_rot == 1) then
     njiyu = 3 * kmp - 6
  else
     njiyu = 3 * kmp
  end if

  temp_now = 2.0e0_PREC * veloet / real(njiyu,PREC)
  rtarget = sqrt(1.0e0_PREC + velo_adjst * (ttemp / temp_now - 1.0e0_PREC))

  do imp = 1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     velo_mp(1:3, imp, irep) = velo_mp(1:3, imp, irep) * rtarget
  end do
  
end subroutine simu_velo_settemp
