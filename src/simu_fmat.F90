! simu_fmat
!> @brief This subroutine is to perform the fluctuation matching.

! See mloop_fmat.F90 for references of the fluctuation matching method.
subroutine simu_fmat()

   use const_physical
   use var_fmat,   only : bl_sum, bl_sum2, ba_sum, ba_sum2,   &
                          dih_sum_A, dih_sum2_A, &
                          dih_sum_B, dih_sum2_B, &
                          nl_sum, nl_sum2, &
                          bp_sum, bp_sum2, &
                          st_sum, st_sum2, &
                          i_num_sum, infmat
   use var_struct, only : xyz_mp_rep,             &
                          nbd, nba, ndih, ncon,   &
                          nrna_bp, nrna_st,       &
                          ibd2mp, iba2mp,         &
                          idih2mp, icon2mp,       &
                          irna_bp2mp, irna_st2mp, &
                          imp2unit, iunit2ba, iunit2dih
   implicit none

   ! local variables
   integer :: idx
   integer :: imp1, imp2, imp3, imp4
   !integer :: iunit1
   real(PREC) :: c11, c12, c13, c21, c22, c23, c33
   real(PREC) :: t1, t3, t4, t3t4
   real(PREC) :: co_dih, co_theta
   real(PREC) :: v21(3), v32(3), v43(3)
   real(PREC) :: dist2, theta, dih
   real(PREC) :: coord_system
   integer, parameter :: IREP = 1

   ! -------------------------------------------------------------------
   i_num_sum = i_num_sum + 1

   ! -------------------------------------------------------------------
   ! bond length
   do idx = 1, nbd
      imp1 = ibd2mp(1, idx)
      imp2 = ibd2mp(2, idx)

      v21(1:3) = xyz_mp_rep(1:3, imp2,IREP) - xyz_mp_rep(1:3, imp1,IREP)
      dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

      bl_sum(idx)  = bl_sum(idx)  + sqrt(dist2)
      bl_sum2(idx) = bl_sum2(idx) + dist2
   enddo

   ! -------------------------------------------------------------------
   ! bond angle
   do idx = 1, nba
      imp1 = iba2mp(1, idx)
      imp2 = iba2mp(2, idx)
      imp3 = iba2mp(3, idx)

      v21(1:3) = xyz_mp_rep(1:3, imp2, IREP) - xyz_mp_rep(1:3, imp1, IREP)
      v32(1:3) = xyz_mp_rep(1:3, imp3, IREP) - xyz_mp_rep(1:3, imp2, IREP)
     
      c11 = v21(1) * v21(1) + v21(2) * v21(2) + v21(3) * v21(3)
      c22 = v32(1) * v32(1) + v32(2) * v32(2) + v32(3) * v32(3)
      c21 = v32(1) * v21(1) + v32(2) * v21(2) + v32(3) * v21(3)
      
      co_theta = - c21 / sqrt(c11 * c22)
     
      if(co_theta > 1.0e0_PREC) then
         co_theta = 1.0e0_PREC
      else if(co_theta < -1.0e0_PREC) then
         co_theta = -1.0e0_PREC
      end if

      theta = acos(co_theta)
      ba_sum(idx)  = ba_sum(idx)  + theta
      ba_sum2(idx) = ba_sum2(idx) + theta*theta
   enddo

   ! -------------------------------------------------------------------
   ! dih angle
   do idx = 1, ndih

      imp1 = idih2mp(1, idx)
      imp2 = idih2mp(2, idx)
      imp3 = idih2mp(3, idx)
      imp4 = idih2mp(4, idx)
 
      v21(1:3) = xyz_mp_rep(1:3, imp2, IREP) - xyz_mp_rep(1:3, imp1, IREP) 
      v32(1:3) = xyz_mp_rep(1:3, imp3, IREP) - xyz_mp_rep(1:3, imp2, IREP)
      v43(1:3) = xyz_mp_rep(1:3, imp4, IREP) - xyz_mp_rep(1:3, imp3, IREP)
      
      c11 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)
      c22 = v32(1)*v32(1) + v32(2)*v32(2) + v32(3)*v32(3)
      c33 = v43(1)*v43(1) + v43(2)*v43(2) + v43(3)*v43(3)
 
      c12 = v21(1)*v32(1) + v21(2)*v32(2) + v21(3)*v32(3)
      c13 = v21(1)*v43(1) + v21(2)*v43(2) + v21(3)*v43(3)
      c23 = v32(1)*v43(1) + v32(2)*v43(2) + v32(3)*v43(3)
 
      t1 = c12*c23 - c13*c22
      t3 = c11*c22 - c12*c12
      t4 = c22*c33 - c23*c23
 
      if (t3 < ZERO_JUDGE) t3 = ZERO_JUDGE
      if (t4 < ZERO_JUDGE) t4 = ZERO_JUDGE
 
      t3t4 = sqrt(t3*t4)
      co_dih      = t1  / t3t4
 
      ! when co_dih > 1 ,or < 1 
      if(co_dih > 1.0e0_PREC) then
         co_dih = 1.0e0_PREC
      else if(co_dih < -1.0e0_PREC) then
         co_dih = -1.0e0_PREC
      end if

      dih = acos(co_dih)
      if (dih > F_PI) dih = F_PI
 
      coord_system = v21(1) * v32(2) * v43(3) + &
                     v32(1) * v43(2) * v21(3) + &
                     v43(1) * v21(2) * v32(3) - &
                     v43(1) * v32(2) * v21(3) - &
                     v21(1) * v43(2) * v32(3) - &
                     v32(1) * v21(2) * v43(3)
 
      !si_dih = sqrt(1.0e0_PREC - co_dih**2)
      !if(si_dih < ZERO_JUDGE) si_dih = ZERO_JUDGE
 
      if(coord_system < 0.0e0_PREC) then
         !si_dih = -si_dih
         dih = -dih
      end if
      
      dih_sum_A(idx)  = dih_sum_A(idx)  + dih
      dih_sum2_A(idx) = dih_sum2_A(idx) + dih*dih

      if (dih < 0.0e0_PREC) then
         dih = dih + 2 * F_PI
      endif
      dih_sum_B(idx)  = dih_sum_B(idx)  + dih
      dih_sum2_B(idx) = dih_sum2_B(idx) + dih*dih
   enddo

   ! -------------------------------------------------------------------
   ! nonlocal
   do idx = 1, ncon
      imp1 = icon2mp(1, idx)
      imp2 = icon2mp(2, idx)
 
      v21(1:3) = xyz_mp_rep(1:3, imp2, IREP) - xyz_mp_rep(1:3, imp1, IREP)
      dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

      nl_sum(idx)  = nl_sum(idx)  + sqrt(dist2)
      nl_sum2(idx) = nl_sum2(idx) + dist2
   enddo

   ! -------------------------------------------------------------------
   ! Base pair
   do idx = 1, nrna_bp
      imp1 = irna_bp2mp(1, idx)
      imp2 = irna_bp2mp(2, idx)
 
      v21(1:3) = xyz_mp_rep(1:3, imp2, IREP) - xyz_mp_rep(1:3, imp1, IREP)
      dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

      bp_sum(idx)  = bp_sum(idx)  + sqrt(dist2)
      bp_sum2(idx) = bp_sum2(idx) + dist2
   enddo

   ! -------------------------------------------------------------------
   ! Base stack
   do idx = 1, nrna_st
      imp1 = irna_st2mp(1, idx)
      imp2 = irna_st2mp(2, idx)
 
      v21(1:3) = xyz_mp_rep(1:3, imp2, IREP) - xyz_mp_rep(1:3, imp1, IREP)
      dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

      st_sum(idx)  = st_sum(idx)  + sqrt(dist2)
      st_sum2(idx) = st_sum2(idx) + dist2
   enddo

endsubroutine simu_fmat
