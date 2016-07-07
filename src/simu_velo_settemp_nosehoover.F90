! ***********************************************************************
subroutine simu_velo_nosehoover(velo_mp, irep, tempk, velo_yojou)

  use const_maxsize
  use const_physical
  use var_setp,   only : insimu, ifix_mp
  use var_struct, only : nmp_real, cmass_mp
  implicit none
      
  ! --------------------------------------------------------------------
  real(PREC), intent(in) :: tempk
  integer,    intent(in) :: irep
  real(PREC), intent(out) :: velo_yojou
  real(PREC), intent(inout) :: velo_mp(:,:,:)

  ! --------------------------------------------------------------------
  ! local variables
  integer :: imp, kmp
  integer :: njiyu
  real(PREC) :: ttemp, veloet2

  ! --------------------------------------------------------------------
  ttemp = tempk * BOLTZ_KCAL_MOL

  veloet2 = 0.0e0_PREC
  kmp = 0
  do imp = 1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     veloet2 = veloet2 + cmass_mp(imp) * (velo_mp(1, imp, irep)**2 &
          + velo_mp(2, imp, irep)**2 + velo_mp(3, imp, irep)**2)
     kmp = kmp + 1
  end do

  
  if(insimu%i_no_trans_rot == 1) then
     njiyu = 3 * kmp - 6
  else
     njiyu = 3 * kmp
  end if

  velo_yojou = veloet2 - njiyu * ttemp
  
end subroutine simu_velo_nosehoover
