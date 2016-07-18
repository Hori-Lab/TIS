! simu_velo_adjst
!> @brief Correct the velocities to remove system translation and rotation

subroutine simu_velo_adjst(velo_mp, irep)

  use const_maxsize
  use var_setp,   only : fix_mp
  use var_struct, only : nmp_real, xyz_mp_rep, cmass_mp

  implicit none

  ! --------------------------------------------------------------------
  real(PREC), intent(inout) :: velo_mp(:,:,:)
  integer,    intent(in)    :: irep

  ! --------------------------------------------------------------------
  integer :: i, j, imp, kmp
  real(PREC) :: cx, cy, cz, hx, hy, hz, sumx, sumy, sumz
  real(PREC) :: txx, txy, txz, tyy, tyz, tzz
  real(PREC) :: anglmt(3), rotate(3, 3), etator(3, 3)
  
  ! --------------------------------------------------------------------
  ! constrain the translation
  ! omega = I^-1 L
  ! L: anglmt
  ! I: rotate
  ! I^-1: etator
  ! cx: omega
  ! --------------------------------------------------------------------

  sumx = 0.0e0_PREC
  sumy = 0.0e0_PREC
  sumz = 0.0e0_PREC

  kmp = 0
  do imp = 1, nmp_real
     if(fix_mp(imp)) cycle
     sumx = velo_mp(1, imp, irep) + sumx
     sumy = velo_mp(2, imp, irep) + sumy
     sumz = velo_mp(3, imp, irep) + sumz
     kmp = kmp + 1
  end do
      
  hx = sumx / real(kmp, PREC)
  hy = sumy / real(kmp, PREC)
  hz = sumz / real(kmp, PREC)

  do imp = 1, nmp_real
     if(fix_mp(imp)) cycle
     velo_mp(1, imp, irep) = velo_mp(1, imp, irep) - hx
     velo_mp(2, imp, irep) = velo_mp(2, imp, irep) - hy
     velo_mp(3, imp, irep) = velo_mp(3, imp, irep) - hz
  end do

  ! --------------------------------------------------------------------
  ! constrain the rotation
  anglmt(1) = 0.0e0_PREC
  anglmt(2) = 0.0e0_PREC
  anglmt(3) = 0.0e0_PREC

  ! calc total angl moment
  do imp = 1, nmp_real
     if(fix_mp(imp)) cycle
     anglmt(1) = anglmt(1) + cmass_mp(imp)                                  &
         * (  xyz_mp_rep(2, imp, irep) * velo_mp(3, imp, irep)      &
            - xyz_mp_rep(3, imp, irep) * velo_mp(2, imp, irep)    )
     anglmt(2) = anglmt(2) + cmass_mp(imp)                                  &
         * (  xyz_mp_rep(3, imp, irep) * velo_mp(1, imp, irep)      &
            - xyz_mp_rep(1, imp, irep) * velo_mp(3, imp, irep)    )
     anglmt(3) = anglmt(3) + cmass_mp(imp)                                  &
         * (  xyz_mp_rep(1, imp, irep) * velo_mp(2, imp, irep)      &
            - xyz_mp_rep(2, imp, irep) * velo_mp(1, imp, irep)    )
  end do
         
  txx = 0.0e0_PREC
  txy = 0.0e0_PREC
  txz = 0.0e0_PREC
  tyy = 0.0e0_PREC
  tyz = 0.0e0_PREC
  tzz = 0.0e0_PREC
      
  do imp = 1, nmp_real
     if(fix_mp(imp)) cycle
     txx = txx + cmass_mp(imp) * xyz_mp_rep(1, imp, irep) * xyz_mp_rep(1, imp, irep)
     txy = txy + cmass_mp(imp) * xyz_mp_rep(1, imp, irep) * xyz_mp_rep(2, imp, irep)
     txz = txz + cmass_mp(imp) * xyz_mp_rep(1, imp, irep) * xyz_mp_rep(3, imp, irep)
     tyy = tyy + cmass_mp(imp) * xyz_mp_rep(2, imp, irep) * xyz_mp_rep(2, imp, irep)
     tyz = tyz + cmass_mp(imp) * xyz_mp_rep(2, imp, irep) * xyz_mp_rep(3, imp, irep)
     tzz = tzz + cmass_mp(imp) * xyz_mp_rep(3, imp, irep) * xyz_mp_rep(3, imp, irep)
  end do
      
  rotate(1, 1) = tyy + tzz
  rotate(2, 1) = -txy
  rotate(3, 1) = -txz
  rotate(1, 2) = -txy
  rotate(2, 2) = txx + tzz
  rotate(3, 2) = -tyz
  rotate(1, 3) = -txz
  rotate(2, 3) = -tyz
  rotate(3, 3) = txx + tyy
       
  ! -------------------------------------------------------------------- 
  do i = 1, 3
     do j = 1, 3
        if(i == j) then
           etator(i, j) = 1.0e0_PREC
        else
           etator(i, j) = 0.0e0_PREC
        end if
     end do
  end do
      
  call util_ppgauss(3, rotate, etator)
         
  cx = etator(1, 1) * anglmt(1) + etator(1, 2) * anglmt(2) + &
       etator(1, 3) * anglmt(3)
  cy = etator(2, 1) * anglmt(1) + etator(2, 2) * anglmt(2) + &
       etator(2, 3) * anglmt(3)
  cz = etator(3, 1) * anglmt(1) + etator(3, 2) * anglmt(2) + &
       etator(3, 3) * anglmt(3)
      
  do imp = 1, nmp_real
     if(fix_mp(imp)) cycle
     velo_mp(1, imp, irep) = velo_mp(1, imp, irep) + &
          cz * xyz_mp_rep(2, imp, irep) - cy * xyz_mp_rep(3, imp, irep)
     velo_mp(2, imp, irep) = velo_mp(2, imp, irep) + &
          cx * xyz_mp_rep(3, imp, irep) - cz * xyz_mp_rep(1, imp, irep)
     velo_mp(3, imp, irep) = velo_mp(3, imp, irep) + &
          cy * xyz_mp_rep(1, imp, irep) - cx * xyz_mp_rep(2, imp, irep)
  enddo

end subroutine simu_velo_adjst
