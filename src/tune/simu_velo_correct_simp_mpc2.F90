! ***********************************************************************
subroutine simu_velo_correct_simp_mpc2(velo_mp, irep)

  use const_maxsize, only : PREC, MXMP, MXSOLV_MPC
  use const_physical, only : BOLTZC
  use var_setp, only : ifix_mp
  use var_struct, only : nmp_real, cmass_mp, xyz_mp_rep
  use var_simu, only : istep, tempk
  use var_mpc, only : inmpc, cmass_solv_mpc, velo_solv_mpc, xyz_solv_mpc
  implicit none
      
  ! --------------------------------------------------------------------
  real(PREC), intent(inout) :: velo_mp(:,:,:)
  integer,    intent(in) :: irep

  ! --------------------------------------------------------------------
  !local variables
  integer :: idimn, is, imp
  integer :: n_solv
  real(PREC) :: P_total_mpc(3), velo_cm_mpc(3), total_mass_mpc
  real(PREC) :: rescale_fact_All, total_kine_E2_All
  integer :: i_count_mp(3), i_count_mp_All
  real(PREC) :: total_moment_interia_ts(3,3)


  ! --------------------------------------------------------------------
  !!initialization
  total_mass_mpc = 0.0e0_PREC
  rescale_fact_All = 0.0e0_PREC
  i_count_mp_All = 0
  total_kine_E2_All = 0.0e0_PREC
  do idimn = 1, 3
     P_total_mpc(idimn) = 0.0e0_PREC
     velo_cm_mpc(idimn) = 0.0e0_PREC
     i_count_mp(idimn) = 0    
  end do
  
  ! --------------------------------------------------------------------
  !!------------------------------------------
  !! total-linear momentum should be zero.
  !!------------------------------------------
  total_mass_mpc = 0.0e0_PREC
  n_solv = inmpc%n_all_solv

  do is = 1, n_solv
     P_total_mpc(1:3) = P_total_mpc(1:3) &
          + cmass_solv_mpc(is) * velo_solv_mpc(1:3, is)
     total_mass_mpc = total_mass_mpc + cmass_solv_mpc(is)
  end do
  
  do imp = 1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     P_total_mpc(1:3) = P_total_mpc(1:3) &
          + cmass_mp(imp) * velo_mp(1:3, imp, irep)
     total_mass_mpc = total_mass_mpc + cmass_mp(imp)
  end do
  
  velo_cm_mpc(1:3) = P_total_mpc(1:3) / total_mass_mpc
  
  !!velocity correction
  n_solv = inmpc%n_all_solv
  do is = 1, n_solv
     velo_solv_mpc(1:3, is) = velo_solv_mpc(1:3, is) - velo_cm_mpc(1:3)
  end do

  do imp = 1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     velo_mp(1:3, imp, irep) = velo_mp(1:3, imp, irep) - velo_cm_mpc(1:3)
  end do
  
     
  if(mod(istep, inmpc%nratio_vcorrect_step) == 0) then

     !!------------------------------------------
     !!calculate rescaling factor for solvent particle velocity to reset T
     !!-----------------------------------------
     i_count_mp_All = 0
     total_kine_E2_All = 0.0e0_PREC 
        
     !!<<for mpc solvent particle>>   
     n_solv = inmpc%n_all_solv
     do is = 1, n_solv
        total_kine_E2_All = total_kine_E2_All + cmass_solv_mpc(is) &
             * (velo_solv_mpc(1, is)**2 &
             + velo_solv_mpc(2, is)**2 + velo_solv_mpc(3, is)**2)
     end do
     i_count_mp_All = i_count_mp_All + 3*n_solv

     rescale_fact_All = sqrt((BOLTZC*tempk)*(i_count_mp_All - 3)/total_kine_E2_All)
     
     !!-----------------------------------------
     !! velocity correct
     !!-----------------------------------------
     n_solv = inmpc%n_all_solv
     do is = 1, n_solv
        velo_solv_mpc(1:3, is) = velo_solv_mpc(1:3, is) * rescale_fact_All
     end do
  end if

end subroutine simu_velo_correct_simp_mpc2
