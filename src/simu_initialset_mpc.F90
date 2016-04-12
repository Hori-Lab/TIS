! ***********************************************************************
! set initial position and velocity for mpc-solvent-particle
! set initial velocity for system-particle (protein,DNA,lipid) 
! ***********************************************************************
subroutine simu_initialset_mpc(tempk_in)

  use const_maxsize, only : PREC, MXMP
  use const_physical, only : BOLTZC
  use const_index
  use var_setp, only : irand, ifix_mp
  use var_struct, only : nmp_real, cmass_mp
  use var_simu, only : velo_mp
  use var_mpc, only : inmpc, cmass_solv_mpc, velo_solv_mpc,xyz_solv_mpc,pbox_origin_mpc,pbox_size_mpc
  use var_replica, only : flg_rep, rep2val, irep2grep, &
                          n_replica_all, n_replica_mpi


  implicit none
  ! --------------------------------------------------------------------
  real(PREC), intent(in) :: tempk_in
  
  ! --------------------------------------------------------------------
  ! function
  !! real(PREC) :: recipe_gasdev
  real(PREC) :: recipe_rand01
  ! --------------------------------------------------------------------
  ! local variables
  integer :: imp, idimn, is, irep, grep
  integer :: n_solv
  integer :: i_free_mpc(3) !! degree of freedom in mpc
  real(PREC) :: P_total(3), total_mass  !! total momentum in periodic boundary box(x,y,z)
  real(PREC) :: v_max, tempk
  real(PREC) :: av_velo(3)
  real(PREC) :: rescale_fact(3)
  real(PREC) :: total_kine_E2(3)

  ! --------------------------------------------------------------------
  ! for the same results
  !do imp = 1, MXNOISE
  !   coef = recipe_gasdev(irand)
  !end do
 
  tempk = tempk_in

  do irep = 1, n_replica_mpi

     grep = irep2grep(irep)

     if (flg_rep(REPTYPE%TEMP)) then
        tempk = rep2val(grep, REPTYPE%TEMP)
     endif 

  !!---------------------------------------------------------------------------
  !! initial position setting for mpc solvent
  n_solv=inmpc%n_all_solv
  do is = 1, n_solv
     do idimn = 1, 3
        xyz_solv_mpc(idimn, is) = pbox_origin_mpc(idimn) + pbox_size_mpc(idimn)*recipe_rand01(irand)
     end do
     !!write (*,*) "inira",xyz_solv_mpc(1, is),xyz_solv_mpc(2, is),xyz_solv_mpc(3, is)
  end do
  !!----------------------------------------------------------------
 
  
  do idimn=1, 3
     P_total(idimn)=0.0e0_PREC
  end do
  total_mass=0.0e0_PREC

  !! Initial velocity setting (for mpc solvent particle)  before correction
  n_solv=inmpc%n_all_solv
  do is = 1, n_solv
     v_max= sqrt(3.0e0_PREC*BOLTZC*tempk/cmass_solv_mpc(is))
     do idimn=1, 3
        velo_solv_mpc(idimn, is)=v_max*(2.0e0_PREC * recipe_rand01(irand) -1.0e0_PREC)
        P_total(idimn)=P_total(idimn)+velo_solv_mpc(idimn, is)*cmass_solv_mpc(is) !! calculation :total momentum
     end do
     total_mass=total_mass+cmass_solv_mpc(is)
  end do

  !! Initial velocity setting (for system particle) before correction
  do imp = 1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     v_max= sqrt(3.0d0*BOLTZC*tempk/cmass_mp(imp))
     do idimn=1, 3
        velo_mp(idimn, imp, irep)=v_max*(2.0e0_PREC * recipe_rand01(irand) -1.0e0_PREC)
        P_total(idimn)=P_total(idimn)+velo_mp(idimn, imp, irep)*cmass_mp(imp) !! calculation :total momentum
     end do
     total_mass=total_mass+cmass_mp(imp)
  end do
  
  !! average velocity
  do idimn = 1, 3
     av_velo(idimn)=P_total(idimn)/total_mass
  end do

  !!********************************************************************
  !! actually angular&linear momentum should be set zero.
  !! However, now skip. 05/12
  !!********************************************************************

  !!*******************************************************************::
  !! calculation rescaling factor for velocity !! total (linear) momentum ==>0, system temperature==>T
  

  do idimn = 1, 3
     i_free_mpc(idimn)=0
     total_kine_E2(idimn)=0.0e0_PREC 
     
     !! <<for mpc solvent particle>>   
     n_solv = inmpc%n_all_solv
     do is = 1, n_solv
        total_kine_E2(idimn) = total_kine_E2(idimn) + (velo_solv_mpc(idimn, is)-av_velo(idimn)) * (velo_solv_mpc(idimn, is)-av_velo(idimn)) * cmass_solv_mpc(is)
        i_free_mpc(idimn)=i_free_mpc(idimn)+1
     end do
     
     !!  <<for system mass particle>>
     do imp = 1, nmp_real
        if(ifix_mp(imp) == 1) cycle
        total_kine_E2(idimn) = total_kine_E2(idimn) + (velo_mp(idimn, imp, irep)-av_velo(idimn)) * (velo_mp(idimn, imp, irep)-av_velo(idimn)) * cmass_mp(imp)
        i_free_mpc(idimn)=i_free_mpc(idimn)+1
     end do
     
     !!rescaling-factor : idimn -axis
     rescale_fact(idimn) = sqrt((BOLTZC*tempk)*(i_free_mpc(idimn)*1.0e0_PREC)/total_kine_E2(idimn))
  end do
  
  
  !!*******************************************************::
  !! velocity correction by rescaling factor, shift velocity
  !! <<for mpc solvent particle>>   
  do is = 1, n_solv
     do idimn=1, 3
        velo_solv_mpc(idimn, is)=rescale_fact(idimn)*(velo_solv_mpc(idimn, is)-av_velo(idimn))
     end do
  end do
  
  !!  <<for system mass particle>>
  do imp = 1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     do idimn=1, 3
        velo_mp(idimn, imp, irep)=rescale_fact(idimn)*(velo_mp(idimn, imp, irep)-av_velo(idimn))
     end do
  end do
  
  
  !!**********************************************************************
  !! check for velocity and energy
  !! check for total momentum
  do idimn=1, 3
     P_total(idimn)=0.0e0_PREC
     total_kine_E2(idimn)=0.0e0_PREC
     i_free_mpc(idimn)=0
  end do
  do imp = 1, nmp_real
     if(ifix_mp(imp) == 1) cycle
     do idimn=1, 3
        i_free_mpc(idimn)=i_free_mpc(idimn)+1
        P_total(idimn)=P_total(idimn)+velo_mp(idimn, imp, irep)*cmass_mp(imp)
        total_kine_E2(idimn)=total_kine_E2(idimn)+velo_mp(idimn, imp, irep)*velo_mp(idimn, imp, irep)*cmass_mp(imp)
     end do
  end do
  n_solv=inmpc%n_all_solv
  do is = 1, n_solv
     do idimn=1, 3
        i_free_mpc(idimn)=i_free_mpc(idimn)+1
        P_total(idimn)=P_total(idimn)+velo_solv_mpc(idimn, is)*cmass_solv_mpc(is)
        !! calculation: total momentum
        total_kine_E2(idimn)=total_kine_E2(idimn)+ velo_solv_mpc(idimn, is)*velo_solv_mpc(idimn, is)*cmass_solv_mpc(is)
     end do
  end do
  write(*,*) "check_result total_P (x,y,z):",P_total(1),P_total(2),P_total(3)
  do idimn=1, 3
     total_kine_E2(idimn)=total_kine_E2(idimn)/(i_free_mpc(idimn)*1.0e0_PREC)
  end do
  write(*,*) "check_result total_E (x,y,z)& kB_T:",total_kine_E2(1),total_kine_E2(2),total_kine_E2(3),",kbt=",tempk * BOLTZC

  end do
  
end subroutine simu_initialset_mpc
