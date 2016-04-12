! ***********************************************************************
subroutine simu_velo_correct_simp_mpc(velo_mp, irep)

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
  real(PREC) :: rescale_fact(3),rescale_fact_All
  real(PREC) :: total_kine_E2(3),total_kine_E2_All
  integer :: i_count_mp(3),i_count_mp_All
  real(PREC) :: AngL_total_mpc(3)
  real(PREC) :: rmx2,rmy2,rmz2,rmxy,rmxz,rmyz !! component of tensor (total moment of interia)
  real(PREC) :: center_xyz_mpc(3)
  real(PREC) :: dis_xyz_solv_mpc(3, MXSOLV_MPC)!,dis_xyz_mp_reduced_mpc(3, MXMP)
  real(PREC) :: dis_xyz_mp_mpc(3, MXMP)

  !!for check
  real(PREC) :: check_P_total(3), check_AL_total(3),check_AL_total2(3)
  real(PREC) :: check_AL_solv(3)
  real(PREC) :: ppkkkk(3)
  real(PREC) :: chch_E
  integer :: chch

  !!initialization
  total_mass_mpc=0.0e0_PREC
  rescale_fact_All=0.0e0_PREC
  i_count_mp_All=0
  total_kine_E2_All=0.0e0_PREC
  do idimn=1, 3
     P_total_mpc(idimn)=0.0e0_PREC
     velo_cm_mpc(idimn)=0.0e0_PREC
     rescale_fact(idimn)=0.0e0_PREC
     total_kine_E2(idimn)=0.0e0_PREC
     i_count_mp(idimn)=0    
     AngL_total_mpc(idimn)=0.0e0_PREC
  end do
  rmx2=0.0e0_PREC
  rmy2=0.0e0_PREC
  rmz2=0.0e0_PREC
  rmxy=0.0e0_PREC
  rmxz=0.0e0_PREC
  rmyz=0.0e0_PREC
  
  
  if(inmpc%i_flag_check_mpc ==1)then
     !!<coordinate: center of mass>
     do idimn=1, 3
        center_xyz_mpc(idimn)=0.0e0_PREC
     end do
     n_solv=inmpc%n_all_solv
     do is=1, n_solv
        do idimn=1, 3
           center_xyz_mpc(idimn) = center_xyz_mpc(idimn) + xyz_solv_mpc(idimn,is) * cmass_solv_mpc(is)
        end do
        total_mass_mpc=total_mass_mpc+cmass_solv_mpc(is)
     end do
     do imp=1, nmp_real
        if(ifix_mp(imp) == 1) cycle
        do idimn=1, 3
           !!center_xyz_mpc(idimn) = center_xyz_mpc(idimn) + xyz_mp_reduced_mpc(idimn,imp) * cmass_mp(imp)
           !! angular momentum should be calculated by real-xyz-mpc, not reduced cordinate.
           center_xyz_mpc(idimn) = center_xyz_mpc(idimn) + xyz_mp_rep(idimn, imp, irep) * cmass_mp(imp)
        end do
        total_mass_mpc=total_mass_mpc+cmass_mp(imp)
     end do
     
     do idimn=1, 3
        center_xyz_mpc(idimn) = center_xyz_mpc(idimn)/total_mass_mpc
     end do
     !! displacement from center of mass
     n_solv=inmpc%n_all_solv
     do is=1, n_solv
        do idimn=1, 3
           dis_xyz_solv_mpc(idimn,is) = xyz_solv_mpc(idimn,is) - center_xyz_mpc(idimn)
        end do
     end do
     do imp=1, nmp_real
        if(ifix_mp(imp) == 1) cycle
        do idimn=1, 3
        !!dis_xyz_mp_reduced_mpc(idimn,imp) = xyz_mp_reduced_mpc(idimn,imp) - center_xyz_mpc(idimn)
           dis_xyz_mp_mpc(idimn,imp) = xyz_mp_rep(idimn, imp, irep) - center_xyz_mpc(idimn)
        end do
     end do
  endif

  
  
  ! --------------------------------------------------------------------
  !!------------------------------------------
  !! total-linear momentum should be zero.
  !!------------------------------------------
  total_mass_mpc=0.0e0_PREC
  n_solv=inmpc%n_all_solv
  do is=1, n_solv
     do idimn=1, 3
        P_total_mpc(idimn)=P_total_mpc(idimn)+cmass_solv_mpc(is) * velo_solv_mpc(idimn, is)
     end do
     total_mass_mpc=total_mass_mpc+cmass_solv_mpc(is)
  end do
  
  do imp=1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     do idimn=1, 3
        P_total_mpc(idimn)=P_total_mpc(idimn)+cmass_mp(imp) * velo_mp(idimn, imp, irep)
     end do
     total_mass_mpc=total_mass_mpc+cmass_mp(imp)
  end do
  
!!!  write(*,*)"rep P_total", P_total_mpc(1),P_total_mpc(2),P_total_mpc(3)

  do idimn=1, 3
     velo_cm_mpc(idimn)=P_total_mpc(idimn)/total_mass_mpc
  end do
  
  !!velocity correction
  n_solv=inmpc%n_all_solv
  do is=1, n_solv
     do idimn=1, 3
        velo_solv_mpc(idimn, is) = velo_solv_mpc(idimn, is) - velo_cm_mpc(idimn)
     end do
  end do
  do imp=1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     do idimn=1, 3
        velo_mp(idimn, imp, irep) = velo_mp(idimn, imp, irep)  - velo_cm_mpc(idimn)
     end do
  end do
  
  
  if(inmpc%i_flag_check_mpc ==1)then
     !!***************************************************
     !!  check total linear and angular momentum
     !!**************************************************
     do idimn=1, 3
        check_P_total(idimn)=0.0e0_PREC
        check_AL_total(idimn)=0.0e0_PREC
        check_AL_total2(idimn)=0.0e0_PREC
     end do
     !<solvent particle>
     n_solv=inmpc%n_all_solv
     do is=1, n_solv
        do idimn=1, 3
           check_P_total(idimn) = check_P_total(idimn) + cmass_solv_mpc(is) * velo_solv_mpc(idimn, is)
        end do
        check_AL_total(1) = check_AL_total(1) + (xyz_solv_mpc(2,is)*velo_solv_mpc(3, is)-xyz_solv_mpc(3,is)*velo_solv_mpc(2, is))*cmass_solv_mpc(is)
        check_AL_total(2) = check_AL_total(2) + (xyz_solv_mpc(3,is)*velo_solv_mpc(1, is)-xyz_solv_mpc(1,is)*velo_solv_mpc(3, is))*cmass_solv_mpc(is)
        check_AL_total(3) = check_AL_total(3) + (xyz_solv_mpc(1,is)*velo_solv_mpc(2, is)-xyz_solv_mpc(2,is)*velo_solv_mpc(1, is))*cmass_solv_mpc(is)
        
        check_AL_total2(1) = check_AL_total2(1) + (dis_xyz_solv_mpc(2,is)*velo_solv_mpc(3, is)-dis_xyz_solv_mpc(3,is)*velo_solv_mpc(2, is))*cmass_solv_mpc(is)
        check_AL_total2(2) = check_AL_total2(2) + (dis_xyz_solv_mpc(3,is)*velo_solv_mpc(1, is)-dis_xyz_solv_mpc(1,is)*velo_solv_mpc(3, is))*cmass_solv_mpc(is)
        check_AL_total2(3) = check_AL_total2(3) + (dis_xyz_solv_mpc(1,is)*velo_solv_mpc(2, is)-dis_xyz_solv_mpc(2,is)*velo_solv_mpc(1, is))*cmass_solv_mpc(is)
     end do
     !!  write(*,*) "bf-solv check_total_AngLxyz",check_AL_total(1),check_AL_total(2),check_AL_total(3)
     !!  write(*,*) "bf-solv check_total22_AngLxyz",check_AL_total2(1),check_AL_total2(2),check_AL_total2(3)
     !!  write(*,*) "bf-solv check_total_Pxyz",check_P_total(1),check_P_total(2),check_P_total(3)
     check_AL_solv(1)=check_AL_total(1)
     check_AL_solv(2)=check_AL_total(2)
     check_AL_solv(3)=check_AL_total(3)
     !<system particle>
     do imp=1, nmp_real
        if(ifix_mp(imp) == 1) cycle
        do idimn=1, 3
           check_P_total(idimn) = check_P_total(idimn) + cmass_mp(imp) * velo_mp(idimn, imp, irep)
        end do
        !!check_AL_total(1) = check_AL_total(1) + (xyz_mp_reduced_mpc(2,imp)*velo_mp(3,imp)-xyz_mp_reduced_mpc(3,imp)*velo_mp(2,imp))*cmass_mp(imp)
        !!check_AL_total(2) = check_AL_total(2) + (xyz_mp_reduced_mpc(3,imp)*velo_mp(1,imp)-xyz_mp_reduced_mpc(1,imp)*velo_mp(3,imp))*cmass_mp(imp)
        !!check_AL_total(3) = check_AL_total(3) + (xyz_mp_reduced_mpc(1,imp)*velo_mp(2,imp)-xyz_mp_reduced_mpc(2,imp)*velo_mp(1,imp))*cmass_mp(imp)
        
        check_AL_total(1) = check_AL_total(1) + (xyz_mp_rep(2,imp,irep)*velo_mp(3,imp,irep)-xyz_mp_rep(3,imp,irep)*velo_mp(2,imp,irep))*cmass_mp(imp)
        check_AL_total(2) = check_AL_total(2) + (xyz_mp_rep(3,imp,irep)*velo_mp(1,imp,irep)-xyz_mp_rep(1,imp,irep)*velo_mp(3,imp,irep))*cmass_mp(imp)
        check_AL_total(3) = check_AL_total(3) + (xyz_mp_rep(1,imp,irep)*velo_mp(2,imp,irep)-xyz_mp_rep(2,imp,irep)*velo_mp(1,imp,irep))*cmass_mp(imp)
        
        !!check_AL_total2(1) = check_AL_total2(1) + (dis_xyz_mp_reduced_mpc(2,imp)*velo_mp(3,imp)-dis_xyz_mp_reduced_mpc(3,imp)*velo_mp(2,imp))*cmass_mp(imp)
        !!check_AL_total2(2) = check_AL_total2(2) + (dis_xyz_mp_reduced_mpc(3,imp)*velo_mp(1,imp)-dis_xyz_mp_reduced_mpc(1,imp)*velo_mp(3,imp))*cmass_mp(imp)
        !!check_AL_total2(3) = check_AL_total2(3) + (dis_xyz_mp_reduced_mpc(1,imp)*velo_mp(2,imp)-dis_xyz_mp_reduced_mpc(2,imp)*velo_mp(1,imp))*cmass_mp(imp)
        
        check_AL_total2(1) = check_AL_total2(1) + (dis_xyz_mp_mpc(2,imp)*velo_mp(3,imp,irep)-dis_xyz_mp_mpc(3,imp)*velo_mp(2,imp,irep))*cmass_mp(imp)
        check_AL_total2(2) = check_AL_total2(2) + (dis_xyz_mp_mpc(3,imp)*velo_mp(1,imp,irep)-dis_xyz_mp_mpc(1,imp)*velo_mp(3,imp,irep))*cmass_mp(imp)
        check_AL_total2(3) = check_AL_total2(3) + (dis_xyz_mp_mpc(1,imp)*velo_mp(2,imp,irep)-dis_xyz_mp_mpc(2,imp)*velo_mp(1,imp,irep))*cmass_mp(imp)
     end do
     !!  write(*,*) "bf-total check_total_AngLxyz",check_AL_total(1),check_AL_total(2),check_AL_total(3)
     !!  write(*,*) "bf-total check_total22_AngLxyz",check_AL_total2(1),check_AL_total2(2),check_AL_total2(3)
     !!  write(*,*) "bf-total check_total_Pxyz",check_P_total(1),check_P_total(2),check_P_total(3)
     !!  write(*,*) "bf-system check_total_AngLxyz",check_AL_total(1)-check_AL_solv(1),check_AL_total(2)-check_AL_solv(2),check_AL_total(3)-check_AL_solv(3)
     !  write(*,*) "bf-total check_total22_AngLxyz",check_AL_total2(1),check_AL_total2(2),check_AL_total2(3)
     !  write(*,*) "bf-total check_total_Pxyz",check_P_total(1),check_P_total(2),check_P_total(3)

     ! --------------------------------------------------------------------
     !!------------------------------------------
     !! total angular-momentum should be zero.
     !!------------------------------------------
     !! deleted
     

     !!***************************************************
     !!  check total linear and angular momentum
     !!**************************************************
     do idimn=1, 3
        check_P_total(idimn)=0.0e0_PREC
        check_AL_total(idimn)=0.0e0_PREC
        check_AL_total2(idimn)=0.0e0_PREC
     end do
     !<solvent particle>
     n_solv=inmpc%n_all_solv
     do is=1, n_solv
        do idimn=1, 3
           check_P_total(idimn) = check_P_total(idimn) + cmass_solv_mpc(is) * velo_solv_mpc(idimn, is)
        end do
        
        check_AL_total(1) = check_AL_total(1) + (dis_xyz_solv_mpc(2,is)*velo_solv_mpc(3, is)-dis_xyz_solv_mpc(3,is)*velo_solv_mpc(2, is))*cmass_solv_mpc(is)
        check_AL_total(2) = check_AL_total(2) + (dis_xyz_solv_mpc(3,is)*velo_solv_mpc(1, is)-dis_xyz_solv_mpc(1,is)*velo_solv_mpc(3, is))*cmass_solv_mpc(is)
        check_AL_total(3) = check_AL_total(3) + (dis_xyz_solv_mpc(1,is)*velo_solv_mpc(2, is)-dis_xyz_solv_mpc(2,is)*velo_solv_mpc(1, is))*cmass_solv_mpc(is)
        
        check_AL_total2(1) = check_AL_total2(1) + (xyz_solv_mpc(2,is)*velo_solv_mpc(3, is)-xyz_solv_mpc(3,is)*velo_solv_mpc(2, is))*cmass_solv_mpc(is)
        check_AL_total2(2) = check_AL_total2(2) + (xyz_solv_mpc(3,is)*velo_solv_mpc(1, is)-xyz_solv_mpc(1,is)*velo_solv_mpc(3, is))*cmass_solv_mpc(is)
        check_AL_total2(3) = check_AL_total2(3) + (xyz_solv_mpc(1,is)*velo_solv_mpc(2, is)-xyz_solv_mpc(2,is)*velo_solv_mpc(1, is))*cmass_solv_mpc(is)
     end do
     
     !!  write(*,*) "af-solv check_total_AngLxyz",check_AL_total(1),check_AL_total(2),check_AL_total(3)
     !!  write(*,*) "af-solv check_total22_AngLxyz",check_AL_total2(1),check_AL_total2(2),check_AL_total2(3)
     ppkkkk(1)=check_AL_total(1)
     ppkkkk(2)=check_AL_total(2)
     ppkkkk(3)=check_AL_total(3)
     
     !<system particle>
     do imp=1, nmp_real
        if(ifix_mp(imp) == 1) cycle
        do idimn=1, 3
           check_P_total(idimn) = check_P_total(idimn) + cmass_mp(imp) * velo_mp(idimn, imp, irep)
        end do
        !!check_AL_total(1) = check_AL_total(1) + (dis_xyz_mp_reduced_mpc(2,imp)*velo_mp(3,imp)-dis_xyz_mp_reduced_mpc(3,imp)*velo_mp(2,imp))*cmass_mp(imp)
        !!check_AL_total(2) = check_AL_total(2) + (dis_xyz_mp_reduced_mpc(3,imp)*velo_mp(1,imp)-dis_xyz_mp_reduced_mpc(1,imp)*velo_mp(3,imp))*cmass_mp(imp)
        !!check_AL_total(3) = check_AL_total(3) + (dis_xyz_mp_reduced_mpc(1,imp)*velo_mp(2,imp)-dis_xyz_mp_reduced_mpc(2,imp)*velo_mp(1,imp))*cmass_mp(imp)
        check_AL_total(1) = check_AL_total(1) + (dis_xyz_mp_mpc(2,imp)*velo_mp(3,imp,irep)-dis_xyz_mp_mpc(3,imp)*velo_mp(2,imp,irep))*cmass_mp(imp)
        check_AL_total(2) = check_AL_total(2) + (dis_xyz_mp_mpc(3,imp)*velo_mp(1,imp,irep)-dis_xyz_mp_mpc(1,imp)*velo_mp(3,imp,irep))*cmass_mp(imp)
        check_AL_total(3) = check_AL_total(3) + (dis_xyz_mp_mpc(1,imp)*velo_mp(2,imp,irep)-dis_xyz_mp_mpc(2,imp)*velo_mp(1,imp,irep))*cmass_mp(imp)
        !!check_AL_total2(1) = check_AL_total2(1) + (xyz_mp_reduced_mpc(2,imp)*velo_mp(3,imp)-xyz_mp_reduced_mpc(3,imp)*velo_mp(2,imp))*cmass_mp(imp)
        !!check_AL_total2(2) = check_AL_total2(2) + (xyz_mp_reduced_mpc(3,imp)*velo_mp(1,imp)-xyz_mp_reduced_mpc(1,imp)*velo_mp(3,imp))*cmass_mp(imp)
        !!check_AL_total2(3) = check_AL_total2(3) + (xyz_mp_reduced_mpc(1,imp)*velo_mp(2,imp)-xyz_mp_reduced_mpc(2,imp)*velo_mp(1,imp))*cmass_mp(imp)
        check_AL_total2(1) = check_AL_total2(1) + (xyz_mp_rep(2,imp,irep)*velo_mp(3,imp,irep)-xyz_mp_rep(3,imp,irep)*velo_mp(2,imp,irep))*cmass_mp(imp)
        check_AL_total2(2) = check_AL_total2(2) + (xyz_mp_rep(3,imp,irep)*velo_mp(1,imp,irep)-xyz_mp_rep(1,imp,irep)*velo_mp(3,imp,irep))*cmass_mp(imp)
        check_AL_total2(3) = check_AL_total2(3) + (xyz_mp_rep(1,imp,irep)*velo_mp(2,imp,irep)-xyz_mp_rep(2,imp,irep)*velo_mp(1,imp,irep))*cmass_mp(imp)
     end do
     !!  write(*,*) "af-total check_total_AngLxyz",check_AL_total(1),check_AL_total(2),check_AL_total(3)
     !!  write(*,*) "af-system check_total_AngLxyz",check_AL_total(1)-ppkkkk(1),check_AL_total(2)-ppkkkk(2),check_AL_total(3)-ppkkkk(3)
     !!  write(*,*) "af check_total_AngLxyz",check_AL_total(1),check_AL_total(2),check_AL_total(3)
     !!  write(*,*) "af check_total22_AngLxyz",check_AL_total2(1),check_AL_total2(2),check_AL_total2(3)
     !!  write(*,*) "af check_total_Pxyz",check_P_total(1),check_P_total(2),check_P_total(3)
  endif
     
if(mod(istep, inmpc%nratio_vcorrect_step) == 0) then	
!!------------------------------------------
!!calculate rescaling factor for solvent particle velocity to reset T
!!-----------------------------------------
  i_count_mp_All=0
  total_kine_E2_All=0.0e0_PREC 
  do idimn = 1, 3
     !!total_kine_E2(idimn)=0.0e0_PREC 
     !! i_count_mp(idimn)=0
     
     !!<<for mpc solvent particle>>   
     n_solv=inmpc%n_all_solv
     do is = 1, n_solv
        total_kine_E2_All = total_kine_E2_All + velo_solv_mpc(idimn, is) * velo_solv_mpc(idimn, is) * cmass_solv_mpc(is)
        !!total_kine_E2(idimn) = total_kine_E2(idimn) + velo_solv_mpc(idimn, is) * velo_solv_mpc(idimn, is) * cmass_solv_mpc(is)
        i_count_mp_All=i_count_mp_All + 1
        !! i_count_mp(idimn) = i_count_mp(idimn) + 1
     end do
     
     !!<<for system mass particle>>
     !!do imp = 1, nmp
     !!   if(ifix_mp(imp) == 1) cycle
     !!   total_kine_E2_All = total_kine_E2_All + velo_mp(idimn, imp) * velo_mp(idimn, imp) * cmass_mp(imp)
     !!   i_count_mp_All=i_count_mp_All + 1
     !!end do
     
     !!i_count_mp_All=i_count_mp_All+i_count_mp(idimn)
     !!total_kine_E2_All=total_kine_E2_All+total_kine_E2(idimn)
     !!rescaling-factor : idimn -axis
     !!rescale_fact(idimn) = sqrt((BOLTZC*tempk)*(i_count_mp(idimn)*1.0e0_PREC)/total_kine_E2(idimn))
  end do
  !!rescale_fact_All=sqrt((BOLTZC*tempk)*((i_count_mp_All - 6) * 1.0e0_PREC)/total_kine_E2_All)
  !!rescale_fact_All=sqrt((BOLTZC*tempk)*((i_count_mp_All - 6) * 1.0e0_PREC)/total_kine_E2_All)
  rescale_fact_All=sqrt((BOLTZC*tempk)*((i_count_mp_All - 3) * 1.0e0_PREC)/total_kine_E2_All)
  !!write (*,*) "i_count_mp_All=",i_count_mp_All
  
!!-----------------------------------------
!! velocity correct
!!-----------------------------------------
  n_solv=inmpc%n_all_solv
  chch=0
  chch_E=0
  do is=1, n_solv
     do idimn=1, 3
        !!velo_solv_mpc(idimn, is) = (velo_solv_mpc(idimn, is) - velo_cm_mpc(idimn)) * rescale_fact_All
        velo_solv_mpc(idimn, is) = velo_solv_mpc(idimn, is) * rescale_fact_All
        !!chch=chch+1
        !!chch_E=chch_E+1.0/2.0*velo_solv_mpc(idimn, is)*velo_solv_mpc(idimn, is)
     end do
  end do
  !!write (*,*) "chcch_E, chch=",chch_E,chch
  !!chch_E=chch_E/chch
  !!write (*,*) "chcch_E=",chch_E

  !!do imp=1, nmp
  !!   if(ifix_mp(imp) == 1) cycle
  !!   do idimn=1, 3
  !!      velo_mp(idimn, imp) = velo_mp(idimn, imp) * rescale_fact_All
  !!   end do
  !!end do

  !!write(*,*) "kana_pp_pass_pp"
  end if
  !!write(*,*) "kana_pp_pass_qq"

end subroutine simu_velo_correct_simp_mpc
