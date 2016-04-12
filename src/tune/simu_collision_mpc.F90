!simu_collision_mpc
!> @brief

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

subroutine simu_collision_mpc(irep)

  use const_maxsize
  use const_physical
!  use const_index
  use if_mloop
!  use if_write
!  use if_energy
  use var_struct, only : nmp_real, xyz_mp_rep, cmass_mp
  use var_simu, only : istep, tstep, velo_mp, tempk
  use var_mpc, only : inmpc, xyz_solv_mpc, velo_solv_mpc, cmass_solv_mpc, &
                      pbox_origin_mpc, pbox_size_mpc

  use time, only : tm_grid_mpc, tm_rotate_mpc, tm_velo_mpc, &
                   time_s, time_e

  implicit none

  ! -----------------------------------------------------------------
  integer, intent(in) :: irep


  ! -----------------------------------------------------------------
  ! local variables
  integer :: imp, idimn
  integer :: is, n_solv, n_mpc_multi
  real(PREC) :: check_velo_bf(MXSOLV_MPC)
  real(PREC) :: check_moment_bf(3), check_moment_af(3)
  real(PREC) :: check_ALmoment_bf(3), check_ALmoment_af(3)
  real(PREC) :: check_kE_bf(3), check_kE_af(3), cm_xyz_check(3)
  integer :: n_out_boundary(3), ibin,idis_velo_solv(3,0:400)
  integer :: ivelo, iii, i_total_v(3)
  real(PREC) :: rdist12, rdist23
  real(PREC) :: ppskkdd(3), ppskkdd2(3)

  ! -----------------------------------------------------------------
  !!call simu_velo_adjst_settemp_mpc(velo_mp, tempk)
  if(mod(istep, inmpc%nratio_colli_step) /= 0) then
     return
  end if


  !!write(*,*) "kanada check"
  !!write(*,*)"uuuuuu", inmpc%nratio_colli_step
  n_solv=inmpc%n_all_solv
     
     
  if(inmpc%i_flag_check_mpc ==1)then
     !!check-------------
     do idimn=1, 3
        check_moment_af(idimn)=0.0e0_PREC
        check_ALmoment_af(idimn)=0.0e0_PREC
     end do
     n_solv=inmpc%n_all_solv
     !! for check
     do is=1, n_solv
        do idimn=1, 3
           check_moment_af(idimn)=check_moment_af(idimn)+velo_solv_mpc(idimn,is)*cmass_solv_mpc(is)
        end do
        !!real(PREC) :: check_ALmoment_af(3),check_ALmoment_af(3)
        check_ALmoment_af(1) = check_ALmoment_af(1) + (xyz_solv_mpc(2,is)*velo_solv_mpc(3, is)-xyz_solv_mpc(3,is)*velo_solv_mpc(2, is))*cmass_solv_mpc(is)
        check_ALmoment_af(2) = check_ALmoment_af(2) + (xyz_solv_mpc(3,is)*velo_solv_mpc(1, is)-xyz_solv_mpc(1,is)*velo_solv_mpc(3, is))*cmass_solv_mpc(is)
        check_ALmoment_af(3) = check_ALmoment_af(3) + (xyz_solv_mpc(1,is)*velo_solv_mpc(2, is)-xyz_solv_mpc(2,is)*velo_solv_mpc(1, is))*cmass_solv_mpc(is)
     end do
     ppskkdd(1)= check_ALmoment_af(1)
     ppskkdd(2)= check_ALmoment_af(2)
     ppskkdd(3)= check_ALmoment_af(3)
     ppskkdd2(1)= check_moment_af(1)
     ppskkdd2(2)= check_moment_af(2)
     ppskkdd2(3)= check_moment_af(3)
     write(*,*) "solv just-before-stream total22_AngLx",check_ALmoment_af(1)
     write(*,*) "solv just-before-stream total22_AngLy",check_ALmoment_af(2)
     write(*,*) "solv just-before-stream total22_AngLz",check_ALmoment_af(3)
     write(*,*) "solv just-before-stream total22_LinearLx",check_moment_af(1)
     write(*,*) "solv just-before-stream total22_LinearLy",check_moment_af(2)
     write(*,*) "solv just-before-stream total22_LinearLz",check_moment_af(3)
     do imp=1, nmp_real
        do idimn=1, 3
           check_moment_af(idimn)=check_moment_af(idimn)+velo_mp(idimn,imp,irep)*cmass_mp(imp)
        end do
        check_ALmoment_af(1) = check_ALmoment_af(1) + (xyz_mp_rep(2,imp,irep)*velo_mp(3,imp,irep)-xyz_mp_rep(3,imp,irep)*velo_mp(2,imp,irep))*cmass_mp(imp)
        check_ALmoment_af(2) = check_ALmoment_af(2) + (xyz_mp_rep(3,imp,irep)*velo_mp(1,imp,irep)-xyz_mp_rep(1,imp,irep)*velo_mp(3,imp,irep))*cmass_mp(imp)
        check_ALmoment_af(3) = check_ALmoment_af(3) + (xyz_mp_rep(1,imp,irep)*velo_mp(2,imp,irep)-xyz_mp_rep(2,imp,irep)*velo_mp(1,imp,irep))*cmass_mp(imp)
     end do
     write(*,*) "total just-before-stream total22_AngLx",check_ALmoment_af(1)
     write(*,*) "total just-before-stream total22_AngLy",check_ALmoment_af(2)
     write(*,*) "total just-before-stream total22_AngLz",check_ALmoment_af(3)
     write(*,*) "total just-before-stream total22_LinearLx",check_moment_af(1)
     write(*,*) "total just-before-stream total22_LinearLy",check_moment_af(2)
     write(*,*) "total just-before-stream total22_LinearLz",check_moment_af(3)
     write(*,*) "system just-before-stream total22_AngLx",check_ALmoment_af(1)-ppskkdd(1)
     write(*,*) "system just-before-stream total22_AngLy",check_ALmoment_af(2)-ppskkdd(2)
     write(*,*) "system just-before-stream total22_AngLz",check_ALmoment_af(3)-ppskkdd(3)
     write(*,*) "system just-before-stream total22_LinearLx",check_moment_af(1)-ppskkdd2(1)
     write(*,*) "system just-before-stream total22_LinearLy",check_moment_af(2)-ppskkdd2(2)
     write(*,*) "system just-before-stream total22_LinearLz",check_moment_af(3)-ppskkdd2(3)
     !!check-------------
  endif
  
  !! av_xyz_coord_solv(1:3)=0.0
  !! av2_xyz_coord_solv(1:3)=0.0
  !!------------------------------------------------(A) start
  !! free-streaming process for mpc solvent particle
  n_solv=inmpc%n_all_solv
  n_mpc_multi=inmpc%nratio_colli_step
  do is=1, n_solv
     !!do idimn=1, 3
     xyz_solv_mpc(1:3, is) = xyz_solv_mpc(1:3, is) + velo_solv_mpc(1:3, is) * (inmpc%nratio_colli_step) * tstep
     !!write(*,*) "t_step",tstep	
     !!end do
     
     !! for checking ======================================
     !!av_xyz_coord_solv(1:3) = av_xyz_coord_solv(1:3) + (xyz_solv_mpc(1:3, is) - xyz_initial_coord_solv(1:3, is))
     !!av2_xyz_coord_solv(1:3) = av2_xyz_coord_solv(1:3) + (xyz_solv_mpc(1:3, is) - xyz_initial_coord_solv(1:3, is))*(xyz_solv_mpc(1:3, is) - xyz_initial_coord_solv(1:3, is))
     
     !!samsamsam=200	
     !!if(is==samsamsam)then
     !!av_xyz_coord_solv_tmp(1:3) = av_xyz_coord_solv(1:3)/samsamsam
     !!av2_xyz_coord_solv_tmp(1:3)= av2_xyz_coord_solv(1:3)/samsamsam
     !!write(*,*) "solv200_xyz_ss_avx",av_xyz_coord_solv_tmp(1)
     !!write(*,*) "solv200_xyz_ss_avy",av_xyz_coord_solv_tmp(2)
     !!write(*,*) "solv200_xyz_ss_avz",av_xyz_coord_solv_tmp(3)
     !!write(*,*) "solv200_xyz_ss_av2x",av2_xyz_coord_solv_tmp(1)
     !!write(*,*) "solv200_xyz_ss_av2y",av2_xyz_coord_solv_tmp(2)
     !!write(*,*) "solv200_xyz_ss_av2z",av2_xyz_coord_solv_tmp(3)
     !! endif	
     
     
     !! for checking ======================================
  end do
  !!av_xyz_coord_solv(1:3) =av_xyz_coord_solv(1:3)/n_solv
  !!av2_xyz_coord_solv(1:3)=av2_xyz_coord_solv(1:3)/n_solv
  
  !!do is=1, 200
  !! for checking ======================================
  !!    av_xyz_coord_solv(1:3) = av_xyz_coord_solv(1:3) + (xyz_solv_mpc(1:3, is) - xyz_initial_coord_solv(1:3, is))
  !!    av2_xyz_coord_solv(1:3) = av2_xyz_coord_solv(1:3) + (xyz_solv_mpc(1:3, is) - xyz_initial_coord_solv(1:3, is))*(xyz_solv_mpc(1:3, is) - xyz_initial_coord_solv(1:3, is))
  !! for checking ======================================
  !!end do	
  !!av_xyz_coord_solv(1:3) =av_xyz_coord_solv(1:3)/200
  !!av2_xyz_coord_solv(1:3)=av2_xyz_coord_solv(1:3)/200
     
  !! write(*,*) "solv_xyz_ss_avx",av_xyz_coord_solv(1)
  !! write(*,*) "solv_xyz_ss_avy",av_xyz_coord_solv(2)
  !! write(*,*) "solv_xyz_ss_avz",av_xyz_coord_solv(3)
  
  !! write(*,*) "solv_xyz_ss_av2x",av2_xyz_coord_solv(1)
  !! write(*,*) "solv_xyz_ss_av2y",av2_xyz_coord_solv(2)
  !! write(*,*) "solv_xyz_ss_av2z",av2_xyz_coord_solv(3)
  
  
  !!write(*,*) "solv_mass", cmass_solv_mpc(1)
  !!write(*,*) "systm_mass", cmass_mp(1)
  !! write(*,*) "solv_xyz_ss_x", xyz_solv_mpc(1, 1)
  !! write(*,*) "solv_xyz_ss_y", xyz_solv_mpc(2, 1)
  !! write(*,*) "solv_xyz_ss_z", xyz_solv_mpc(3, 1)
  
  !! write(*,*) "system_xyz_x",  xyz_mp(1, 1)
  !! write(*,*) "system_xyz_y",  xyz_mp(2, 1)
  !! write(*,*) "system_xyz_z",  xyz_mp(3, 1)
  
  !!-----------------------------------------------(A) end
  
  
  if(inmpc%i_flag_check_mpc ==1)then
     !!check-------------
     do idimn=1, 3
        check_moment_af(idimn)=0.0e0_PREC
        check_ALmoment_af(idimn)=0.0e0_PREC
     end do
     n_solv=inmpc%n_all_solv
     !! for check
     do is=1, n_solv
        do idimn=1, 3
           check_moment_af(idimn)=check_moment_af(idimn)+velo_solv_mpc(idimn,is)*cmass_solv_mpc(is)
        end do
        !!real(PREC) :: check_ALmoment_af(3),check_ALmoment_af(3)
        check_ALmoment_af(1) = check_ALmoment_af(1) + (xyz_solv_mpc(2,is)*velo_solv_mpc(3, is)-xyz_solv_mpc(3,is)*velo_solv_mpc(2, is))*cmass_solv_mpc(is)
        check_ALmoment_af(2) = check_ALmoment_af(2) + (xyz_solv_mpc(3,is)*velo_solv_mpc(1, is)-xyz_solv_mpc(1,is)*velo_solv_mpc(3, is))*cmass_solv_mpc(is)
        check_ALmoment_af(3) = check_ALmoment_af(3) + (xyz_solv_mpc(1,is)*velo_solv_mpc(2, is)-xyz_solv_mpc(2,is)*velo_solv_mpc(1, is))*cmass_solv_mpc(is)
     end do
     ppskkdd(1)= check_ALmoment_af(1)
     ppskkdd(2)= check_ALmoment_af(2)
     ppskkdd(3)= check_ALmoment_af(3)
     ppskkdd2(1)= check_moment_af(1)
     ppskkdd2(2)= check_moment_af(2)
     ppskkdd2(3)= check_moment_af(3)
     write(*,*) "solv just-after-stream total22_AngLx",check_ALmoment_af(1)
     write(*,*) "solv just-after-stream total22_AngLy",check_ALmoment_af(2)
     write(*,*) "solv just-after-stream total22_AngLz",check_ALmoment_af(3)
     write(*,*) "solv just-after-stream total22_LinearLx",check_moment_af(1)
     write(*,*) "solv just-after-stream total22_LinearLy",check_moment_af(2)
     write(*,*) "solv just-after-stream total22_LinearLz",check_moment_af(3)
     do imp=1, nmp_real
        do idimn=1, 3
           check_moment_af(idimn)=check_moment_af(idimn)+velo_mp(idimn,imp,irep)*cmass_mp(imp)
        end do
        check_ALmoment_af(1) = check_ALmoment_af(1) + (xyz_mp_rep(2,imp,irep)*velo_mp(3,imp,irep)-xyz_mp_rep(3,imp,irep)*velo_mp(2,imp,irep))*cmass_mp(imp)
        check_ALmoment_af(2) = check_ALmoment_af(2) + (xyz_mp_rep(3,imp,irep)*velo_mp(1,imp,irep)-xyz_mp_rep(1,imp,irep)*velo_mp(3,imp,irep))*cmass_mp(imp)
        check_ALmoment_af(3) = check_ALmoment_af(3) + (xyz_mp_rep(1,imp,irep)*velo_mp(2,imp,irep)-xyz_mp_rep(2,imp,irep)*velo_mp(1,imp,irep))*cmass_mp(imp)
     end do
     write(*,*) "total just-after-stream total22_AngLx",check_ALmoment_af(1)
     write(*,*) "total just-after-stream total22_AngLy",check_ALmoment_af(2)
     write(*,*) "total just-after-stream total22_AngLz",check_ALmoment_af(3)
     write(*,*) "total just-after-stream total22_LinearLx",check_moment_af(1)
     write(*,*) "total just-after-stream total22_LinearLy",check_moment_af(2)
     write(*,*) "total just-after-stream total22_LinearLz",check_moment_af(3)
     write(*,*) "system just-after-stream total22_AngLx",check_ALmoment_af(1)-ppskkdd(1)
     write(*,*) "system just-after-stream total22_AngLy",check_ALmoment_af(2)-ppskkdd(2)
     write(*,*) "system just-after-stream total22_AngLz",check_ALmoment_af(3)-ppskkdd(3)
     write(*,*) "system just-after-stream total22_LinearLx",check_moment_af(1)-ppskkdd2(1)
     write(*,*) "system just-after-stream total22_LinearLy",check_moment_af(2)-ppskkdd2(2)
     write(*,*) "system just-after-stream total22_LinearLz",check_moment_af(3)-ppskkdd2(3)
     !!check-------------
  endif
  
  

  !! for check
  !!do is=1, n_solv
  !!   write(*,*) "xyz_solv_mpc:BF",xyz_solv_mpc(1, is),xyz_solv_mpc(2, is),xyz_solv_mpc(3, is)
  !!end do
  
     
  !!-----------------------------------------------(B) start
  !! grid (cell) devision process
  !! system particle and solvent particle of mpc is devided into cell
  !!--------------
  !! input: xyz_solv_mpc(1:3, MXSOLV_MPC), xyz_mp(1:3, MXMP) :coordinate
  !! output: isolv2grid_mpc(MXSOLV_MPC),imp2grid_mpc(MXMP) : relation between mass_point_id and grid-number
  TIME_S( tm_grid_mpc )
  call simu_grid_devision_mpc(irep)
  TIME_E( tm_grid_mpc )
  !!-----------------------------------------------(B) end
  
  
  !!check-------------
  !!check_ALmoment_af(1)=0.0e0_PREC
  !!check_ALmoment_af(2)=0.0e0_PREC
  !!check_ALmoment_af(3)=0.0e0_PREC
  !!n_solv=inmpc%n_all_solv
  !!do is=1, n_solv
  !!   check_ALmoment_af(1) = check_ALmoment_af(1) + (xyz_solv_mpc(2,is)*velo_solv_mpc(3, is)-xyz_solv_mpc(3,is)*velo_solv_mpc(2, is))*cmass_solv_mpc(is)
  !!   check_ALmoment_af(2) = check_ALmoment_af(2) + (xyz_solv_mpc(3,is)*velo_solv_mpc(1, is)-xyz_solv_mpc(1,is)*velo_solv_mpc(3, is))*cmass_solv_mpc(is)
  !!   check_ALmoment_af(3) = check_ALmoment_bf(3) + (xyz_solv_mpc(1,is)*velo_solv_mpc(2, is)-xyz_solv_mpc(2,is)*velo_solv_mpc(1, is))*cmass_solv_mpc(is)
  !!end do
  !!write(*,*) "solv just-after-grid total22_AngLx",check_ALmoment_af(1)
  !!write(*,*) "solv just-after-grid total22_AngLy",check_ALmoment_af(2)
  !!write(*,*) "solv just-after-grid total22_AngLz",check_ALmoment_af(3)
  !!check-------------
  !! for check
  !!do is=1, n_solv
  !!   write(*,*) "xyz_solv_mpc:AF",xyz_solv_mpc(1, is),xyz_solv_mpc(2, is),xyz_solv_mpc(3, is),isolv2grid_mpc(is)
  !!end do
  
  
  if(inmpc%i_flag_check_mpc ==1)then
     !! for check
     do idimn=1, 3
        check_moment_bf(idimn)=0.0e0_PREC
        check_moment_af(idimn)=0.0e0_PREC
        check_ALmoment_bf(idimn)=0.0e0_PREC
        check_ALmoment_af(idimn)=0.0e0_PREC
        check_kE_bf(idimn)=0.0e0_PREC
        check_kE_af(idimn)=0.0e0_PREC
        n_out_boundary(idimn)=0
        cm_xyz_check(idimn)=0
     end do
     
     !! for check
     do is=1, n_solv
        check_velo_bf(is) = velo_solv_mpc(1,is)*velo_solv_mpc(1,is)+velo_solv_mpc(2,is)*velo_solv_mpc(2,is)+velo_solv_mpc(3,is)*velo_solv_mpc(3,is)
        do idimn=1, 3
           check_moment_bf(idimn)=check_moment_bf(idimn)+velo_solv_mpc(idimn,is)*cmass_solv_mpc(is)
           check_kE_bf(idimn)=check_kE_bf(idimn)+1.0e0_PREC/2.0e0_PREC*cmass_solv_mpc(is)*velo_solv_mpc(idimn,is)*velo_solv_mpc(idimn,is)
        end do
        !!real(PREC) :: check_ALmoment_bf(3),check_ALmoment_af(3)
        check_ALmoment_bf(1) = check_ALmoment_bf(1) + (xyz_solv_mpc(2,is)*velo_solv_mpc(3, is)-xyz_solv_mpc(3,is)*velo_solv_mpc(2, is))*cmass_solv_mpc(is)
        check_ALmoment_bf(2) = check_ALmoment_bf(2) + (xyz_solv_mpc(3,is)*velo_solv_mpc(1, is)-xyz_solv_mpc(1,is)*velo_solv_mpc(3, is))*cmass_solv_mpc(is)
        check_ALmoment_bf(3) = check_ALmoment_bf(3) + (xyz_solv_mpc(1,is)*velo_solv_mpc(2, is)-xyz_solv_mpc(2,is)*velo_solv_mpc(1, is))*cmass_solv_mpc(is)
     end do
     ppskkdd(1)= check_ALmoment_bf(1)
     ppskkdd(2)= check_ALmoment_bf(2)
     ppskkdd(3)= check_ALmoment_bf(3)
     
     ppskkdd2(1)= check_moment_bf(1)
     ppskkdd2(2)= check_moment_bf(2)
     ppskkdd2(3)= check_moment_bf(3)
     
     write(*,*) "solv after-grid total22_AngLx",check_ALmoment_bf(1)
     write(*,*) "solv after-grid total22_AngLy",check_ALmoment_bf(2)
     write(*,*) "solv after-grid total22_AngLz",check_ALmoment_bf(3)
     write(*,*) "solv after-grid total22_LinearLx",check_moment_bf(1)
     write(*,*) "solv after-grid total22_LinearLy",check_moment_bf(2)
     write(*,*) "solv after-grid total22_LinearLz",check_moment_bf(3)
     
     !! for check
     do imp=1, nmp_real
        do idimn=1, 3
           if((xyz_mp_rep(idimn,imp,irep) < pbox_origin_mpc(idimn)) .or.(xyz_mp_rep(idimn,imp,irep) > pbox_origin_mpc(idimn)+pbox_size_mpc(idimn)))then
              n_out_boundary(idimn)=n_out_boundary(idimn)+1
           end if
           cm_xyz_check(idimn)=cm_xyz_check(idimn)+xyz_mp_rep(idimn,imp,irep)
           check_moment_bf(idimn)=check_moment_bf(idimn)+velo_mp(idimn,imp,irep)*cmass_mp(imp)
           check_kE_bf(idimn)=check_kE_bf(idimn)+1.0e0_PREC/2.0e0_PREC*cmass_mp(imp)*velo_mp(idimn,imp,irep)*velo_mp(idimn,imp,irep)
        end do
        check_ALmoment_bf(1) = check_ALmoment_bf(1) + (xyz_mp_rep(2,imp,irep)*velo_mp(3,imp,irep)-xyz_mp_rep(3,imp,irep)*velo_mp(2,imp,irep))*cmass_mp(imp)
        check_ALmoment_bf(2) = check_ALmoment_bf(2) + (xyz_mp_rep(3,imp,irep)*velo_mp(1,imp,irep)-xyz_mp_rep(1,imp,irep)*velo_mp(3,imp,irep))*cmass_mp(imp)
        check_ALmoment_bf(3) = check_ALmoment_bf(3) + (xyz_mp_rep(1,imp,irep)*velo_mp(2,imp,irep)-xyz_mp_rep(2,imp,irep)*velo_mp(1,imp,irep))*cmass_mp(imp)
     end do
     
     write(*,*) "total after-grid total22_AngLx",check_ALmoment_bf(1)
     write(*,*) "total after-grid total22_AngLy",check_ALmoment_bf(2)
     write(*,*) "total after-grid total22_AngLz",check_ALmoment_bf(3)
     write(*,*) "total after-grid total22_LinearLx",check_moment_bf(1)
     write(*,*) "total after-grid total22_LinearLy",check_moment_bf(2)
     write(*,*) "total after-grid total22_LinearLz",check_moment_bf(3)
     
     write(*,*) "system after-grid total22_AngLx",check_ALmoment_bf(1)-ppskkdd(1)
     write(*,*) "system after-grid total22_AngLy",check_ALmoment_bf(2)-ppskkdd(2)
     write(*,*) "system after-grid total22_AngLz",check_ALmoment_bf(3)-ppskkdd(3)
     write(*,*) "system after-grid total22_LinearLx",check_moment_bf(1)-ppskkdd2(1)
     write(*,*) "system after-grid total22_LinearLy",check_moment_bf(2)-ppskkdd2(2)
     write(*,*) "system after-grid total22_LinearLz",check_moment_bf(3)-ppskkdd2(3)
     !!write(*,*) "bf-total before-rotate check_total22_AngLxyz",check_ALmoment_bf(1),check_ALmoment_bf(2),check_ALmoment_bf(3)
     !!write(*,*) "bf-system before-rotate check_total22_AngLxyz",check_ALmoment_bf(1)-ppskkdd(1),check_ALmoment_bf(2)-ppskkdd(2),check_ALmoment_bf(3)-ppskkdd(3)
  endif
  
  !!,velo_mp(1,2),velo_mp(2,2),velo_mp(3,2),velo_mp(1,3),velo_mp(2,3),velo_mp(3,3)
  
  !!--------------------------------------------(C) start
  !! rotate velocity
  TIME_S( tm_rotate_mpc )
!  call simu_rotate_velo_mpc(velo_mp, irep)
  call simu_rotate_velo_mpc_rev(velo_mp, irep)
  TIME_E( tm_rotate_mpc )
  !!--------------------------------------------(C) end
  
  
  if(inmpc%i_flag_check_mpc ==1)then
     !!check
     do idimn=1, 3
        check_moment_bf(idimn)=0.0e0_PREC
        check_moment_af(idimn)=0.0e0_PREC
        check_ALmoment_bf(idimn)=0.0e0_PREC
     end do
     !! for check
     do is=1, n_solv
        check_velo_bf(is) = velo_solv_mpc(1,is)*velo_solv_mpc(1,is)+velo_solv_mpc(2,is)*velo_solv_mpc(2,is)+velo_solv_mpc(3,is)*velo_solv_mpc(3,is)
        do idimn=1, 3
           check_moment_bf(idimn)=check_moment_bf(idimn)+velo_solv_mpc(idimn,is)*cmass_solv_mpc(is)
           !!check_kE_bf(idimn)=check_kE_bf(idimn)+1.0e0_PREC/2.0e0_PREC*cmass_solv_mpc(is)*velo_solv_mpc(idimn,is)*velo_solv_mpc(idimn,is)
        end do
        !!real(PREC) :: check_ALmoment_bf(3),check_ALmoment_af(3)
        check_ALmoment_bf(1) = check_ALmoment_bf(1) + (xyz_solv_mpc(2,is)*velo_solv_mpc(3, is)-xyz_solv_mpc(3,is)*velo_solv_mpc(2, is))*cmass_solv_mpc(is)
        check_ALmoment_bf(2) = check_ALmoment_bf(2) + (xyz_solv_mpc(3,is)*velo_solv_mpc(1, is)-xyz_solv_mpc(1,is)*velo_solv_mpc(3, is))*cmass_solv_mpc(is)
        check_ALmoment_bf(3) = check_ALmoment_bf(3) + (xyz_solv_mpc(1,is)*velo_solv_mpc(2, is)-xyz_solv_mpc(2,is)*velo_solv_mpc(1, is))*cmass_solv_mpc(is)
     end do
     ppskkdd(1)= check_ALmoment_bf(1)
     ppskkdd(2)= check_ALmoment_bf(2)
     ppskkdd(3)= check_ALmoment_bf(3)
     ppskkdd2(1)= check_moment_bf(1)
     ppskkdd2(2)= check_moment_bf(2)
     ppskkdd2(3)= check_moment_bf(3)
     
     write(*,*) "solv after-rotate total22_AngLx",check_ALmoment_bf(1)
     write(*,*) "solv after-rotate total22_AngLy",check_ALmoment_bf(2)
     write(*,*) "solv after-rotate total22_AngLz",check_ALmoment_bf(3)
     write(*,*) "solv after-rotate total22_LinearLx",check_moment_bf(1)
     write(*,*) "solv after-rotate total22_LinearLy",check_moment_bf(2)
     write(*,*) "solv after-rotate total22_LinearLz",check_moment_bf(3)
     
     !! for check
     do imp=1, nmp_real
        do idimn=1, 3
           check_moment_bf(idimn)=check_moment_bf(idimn)+velo_mp(idimn,imp,irep)*cmass_mp(imp)
        end do
        check_ALmoment_bf(1) = check_ALmoment_bf(1) + (xyz_mp_rep(2,imp,irep)*velo_mp(3,imp,irep)-xyz_mp_rep(3,imp,irep)*velo_mp(2,imp,irep))*cmass_mp(imp)
        check_ALmoment_bf(2) = check_ALmoment_bf(2) + (xyz_mp_rep(3,imp,irep)*velo_mp(1,imp,irep)-xyz_mp_rep(1,imp,irep)*velo_mp(3,imp,irep))*cmass_mp(imp)
        check_ALmoment_bf(3) = check_ALmoment_bf(3) + (xyz_mp_rep(1,imp,irep)*velo_mp(2,imp,irep)-xyz_mp_rep(2,imp,irep)*velo_mp(1,imp,irep))*cmass_mp(imp)
     end do
     
     write(*,*) "total after-rotate total22_AngLx",check_ALmoment_bf(1)
     write(*,*) "total after-rotate total22_AngLy",check_ALmoment_bf(2)
     write(*,*) "total after-rotate total22_AngLz",check_ALmoment_bf(3)
     write(*,*) "total after-rotate total22_LinearLx",check_moment_bf(1)
     write(*,*) "total after-rotate total22_LinearLy",check_moment_bf(2)
     write(*,*) "total after-rotate total22_LinearLz",check_moment_bf(3)
     
     write(*,*) "system after-rotate total22_AngLx",check_ALmoment_bf(1)-ppskkdd(1)
     write(*,*) "system after-rotate total22_AngLy",check_ALmoment_bf(2)-ppskkdd(2)
     write(*,*) "system after-rotate total22_AngLz",check_ALmoment_bf(3)-ppskkdd(3)
     write(*,*) "system after-rotate total22_LinearLx",check_moment_bf(1)-ppskkdd2(1)
     write(*,*) "system after-rotate total22_LinearLy",check_moment_bf(2)-ppskkdd2(2)
     write(*,*) "system after-rotate total22_LinearLz",check_moment_bf(3)-ppskkdd2(3)
     !!write(*,*) "bf-total before-rotate check_total22_AngLxyz",check_ALmoment_bf(1),check_ALmoment_bf(2),check_ALmoment_bf(3)
     !!write(*,*) "bf-system before-rotate check_total22_AngLxyz",check_ALmoment_bf(1)-ppskkdd(1),check_ALmoment_bf(2)-ppskkdd(2),check_ALmoment_bf(3)-ppskkdd(3)
  endif
  
  
  !!call simu_velo_correct_settemp_mpc(velo_mp, tempk)
  TIME_S( tm_velo_mpc )
  call simu_velo_correct_simp_mpc(velo_mp, irep)
  TIME_E( tm_velo_mpc )
  
     
  if(inmpc%i_flag_check2_mpc ==1)then
     write(*,*) "kkk1",velo_mp(1,1,irep),velo_mp(2,1,irep),velo_mp(3,1,irep)
     write(*,*) "kkk2",velo_mp(1,2,irep),velo_mp(2,2,irep),velo_mp(3,2,irep)
     write(*,*) "kkk3",velo_mp(1,3,irep),velo_mp(2,3,irep),velo_mp(3,3,irep)
     rdist12 = sqrt((xyz_mp_rep(1, 1, irep)-xyz_mp_rep(1, 2, irep))**2 &
     + (xyz_mp_rep(2, 1, irep)-xyz_mp_rep(2, 2, irep))**2 &
     + (xyz_mp_rep(3, 1, irep)-xyz_mp_rep(3, 2, irep))**2)
     rdist23 = sqrt((xyz_mp_rep(1, 3, irep)-xyz_mp_rep(1, 2, irep))**2 &
     + (xyz_mp_rep(2, 3, irep)-xyz_mp_rep(2, 2, irep))**2 &
     + (xyz_mp_rep(3, 3, irep)-xyz_mp_rep(3, 2, irep))**2)
     write(*,*) "rrr12",rdist12
     write(*,*) "rrr23",rdist23
     write(*,*) "solv001",velo_solv_mpc(1,1),velo_solv_mpc(2,1),velo_solv_mpc(3,1)
     !! linear momentum should be reset (accumulated numerical error in streaming phase).
     !! solvent kinetic energy rescaling for solvent temperature-rescaling T. 
  endif
  
     
     
  if(inmpc%i_flag_check_mpc ==1)then
     !! for check
     do is=1, n_solv
        check_velo_bf(is) = velo_solv_mpc(1,is)*velo_solv_mpc(1,is)+velo_solv_mpc(2,is)*velo_solv_mpc(2,is)+velo_solv_mpc(3,is)*velo_solv_mpc(3,is)
        do idimn=1, 3
           check_moment_af(idimn)=check_moment_af(idimn)+velo_solv_mpc(idimn,is)*cmass_solv_mpc(is)
           check_kE_af(idimn)=check_kE_af(idimn)+1.0e0_PREC/2.0e0_PREC*cmass_solv_mpc(is)*velo_solv_mpc(idimn,is)*velo_solv_mpc(idimn,is)
        end do
        !!real(PREC) :: check_ALmoment_af(3),check_ALmoment_af(3)
        check_ALmoment_af(1) = check_ALmoment_af(1) + (xyz_solv_mpc(2,is)*velo_solv_mpc(3, is)-xyz_solv_mpc(3,is)*velo_solv_mpc(2, is))*cmass_solv_mpc(is)
        check_ALmoment_af(2) = check_ALmoment_af(2) + (xyz_solv_mpc(3,is)*velo_solv_mpc(1, is)-xyz_solv_mpc(1,is)*velo_solv_mpc(3, is))*cmass_solv_mpc(is)
        check_ALmoment_af(3) = check_ALmoment_bf(3) + (xyz_solv_mpc(1,is)*velo_solv_mpc(2, is)-xyz_solv_mpc(2,is)*velo_solv_mpc(1, is))*cmass_solv_mpc(is)
     end do
     !!ppskkdd(1)= check_ALmoment_af(1)
     !!ppskkdd(2)= check_ALmoment_af(2)
     !!ppskkdd(3)= check_ALmoment_af(3)
     !!write(*,*) "af-solv before-rotate check_total22_AngLxyz",check_ALmoment_af(1),check_ALmoment_af(2),check_ALmoment_af(3)
     ppskkdd(1)= check_ALmoment_af(1)
     ppskkdd(2)= check_ALmoment_af(2)
     ppskkdd(3)= check_ALmoment_af(3)
     
     ppskkdd2(1)= check_moment_af(1)
     ppskkdd2(2)= check_moment_af(2)
     ppskkdd2(3)= check_moment_af(3)
     
     write(*,*) "solv after-corect total22_AngLx",check_ALmoment_af(1)
     write(*,*) "solv after-corect total22_AngLy",check_ALmoment_af(2)
     write(*,*) "solv after-corect total22_AngLz",check_ALmoment_af(3)
     write(*,*) "solv after-corect total22_LinearLx",check_moment_af(1)
     write(*,*) "solv after-corect total22_LinearLy",check_moment_af(2)
     write(*,*) "solv after-corect total22_LinearLz",check_moment_af(3)
     
     !! for check
     do imp=1, nmp_real
        !!??????check_velo_bf(is) = velo_solv_mpc(1,is)*velo_solv_mpc(1,is)+velo_solv_mpc(2,is)*velo_solv_mpc(2,is)+velo_solv_mpc(3,is)*velo_solv_mpc(3,is)
        do idimn=1, 3
           check_moment_af(idimn)=check_moment_af(idimn)+velo_mp(idimn,imp,irep)*cmass_mp(imp)
           check_kE_af(idimn)=check_kE_af(idimn)+1.0e0_PREC/2.0e0_PREC*cmass_mp(imp)*velo_mp(idimn,imp,irep)*velo_mp(idimn,imp,irep)
        end do
        check_ALmoment_af(1) = check_ALmoment_af(1) + (xyz_mp_rep(2,imp,irep)*velo_mp(3,imp,irep)-xyz_mp_rep(3,imp,irep)*velo_mp(2,imp,irep))*cmass_mp(imp)
        check_ALmoment_af(2) = check_ALmoment_af(2) + (xyz_mp_rep(3,imp,irep)*velo_mp(1,imp,irep)-xyz_mp_rep(1,imp,irep)*velo_mp(3,imp,irep))*cmass_mp(imp)
        check_ALmoment_af(3) = check_ALmoment_af(3) + (xyz_mp_rep(1,imp,irep)*velo_mp(2,imp,irep)-xyz_mp_rep(2,imp,irep)*velo_mp(1,imp,irep))*cmass_mp(imp)
     end do
     !!write(*,*) "af-total before-rotate check_total22_AngLxyz",check_ALmoment_af(1),check_ALmoment_af(2),check_ALmoment_af(3)
     !!write(*,*) "af-system before-rotate check_total22_AngLxyz",check_ALmoment_af(1)-ppskkdd(1),check_ALmoment_af(2)-ppskkdd(2),check_ALmoment_af(3)-ppskkdd(3)
     
     write(*,*) "total after-corect total22_AngLx",check_ALmoment_af(1)
     write(*,*) "total after-corect total22_AngLy",check_ALmoment_af(2)
     write(*,*) "total after-corect total22_AngLz",check_ALmoment_af(3)
     write(*,*) "total after-corect total22_LinearLx",check_moment_af(1)
     write(*,*) "total after-corect total22_LinearLy",check_moment_af(2)
     write(*,*) "total after-corect total22_LinearLz",check_moment_af(3)
     write(*,*) "system after-corect total22_AngLx",check_ALmoment_af(1)-ppskkdd(1)
     write(*,*) "system after-corect total22_AngLy",check_ALmoment_af(2)-ppskkdd(2)
     write(*,*) "system after-corect total22_AngLz",check_ALmoment_af(3)-ppskkdd(3)
     write(*,*) "system after-corect total22_LinearLx",check_moment_af(1)-ppskkdd2(1)
     write(*,*) "system after-corect total22_LinearLy",check_moment_af(2)-ppskkdd2(2)
     write(*,*) "system after-corect total22_LinearLz",check_moment_af(3)-ppskkdd2(3)
     !! for check
     cm_xyz_check(1)=cm_xyz_check(1)/nmp_real
     cm_xyz_check(2)=cm_xyz_check(2)/nmp_real
     cm_xyz_check(3)=cm_xyz_check(3)/nmp_real
     
     !!write(*,*) "check moment x",check_moment_bf(1),check_moment_af(1)
     !!write(*,*) "check moment y",check_moment_bf(2),check_moment_af(2)
     !!write(*,*) "check moment z",check_moment_bf(3),check_moment_af(3)
     !!write(*,*) "check Kenrgy x",check_kE_bf(1), check_kE_af(1)
     !!write(*,*) "check Kenrgy y",check_kE_bf(2), check_kE_af(2)
     !!write(*,*) "check Kenrgy z",check_kE_bf(3), check_kE_af(3)
     !!write(*,*) "check Kenrgy a",check_kE_bf(1)+check_kE_bf(2)+check_kE_bf(3), check_kE_af(1)+check_kE_af(2)+check_kE_af(3)
     write(*,*) "check Kenegy T",(check_kE_bf(1)+check_kE_bf(2)+check_kE_bf(3))/(nmp_real+n_solv),(check_kE_af(1)+check_kE_af(2)+check_kE_af(3))/(nmp_real+n_solv),&
          3.0e0_PREC/2.0e0_PREC*tempk*BOLTZC
     
  endif
  
  !!write(*,*) "check Outbox xyz",n_out_boundary(1),n_out_boundary(2),n_out_boundary(3)
  !!write(*,*) "check averag xyz",cm_xyz_check(1),cm_xyz_check(2),cm_xyz_check(3)
  
  !!do is=1, n_solv
  !!   check_velo_af(is) = velo_solv_mpc(1,is)*velo_solv_mpc(1,is)+velo_solv_mpc(2,is)*velo_solv_mpc(2,is)+velo_solv_mpc(3,is)*velo_solv_mpc(3,is)
  !!   write(*,*) "BF_solv_velo", check_velo_bf(is),check_velo_af(is)
  !!end do
  
  !!--------------------------------------
  !! momentum, temperature (total-E) collection
  !!-------------------------------------
  !!call simu_velo_adjst_settemp_mpc(velo_mp, tempk)
  
                 
  if(inmpc%i_flag_check2_mpc ==1)then
     if((istep== (inmpc%nratio_colli_step)*1) .or. &
          (istep== (inmpc%nratio_colli_step)*10).or. &
          (istep== (inmpc%nratio_colli_step)*300).or. &
          (istep== (inmpc%nratio_colli_step)*200).or. &
          (istep== (inmpc%nratio_colli_step)*100).or. &
          (istep== (inmpc%nratio_colli_step)*10000).or. &
          (istep== (inmpc%nratio_colli_step)*20000).or. &
          (istep== (inmpc%nratio_colli_step)*30000)) then
        !!check the distribution of velocity
        !! for check
        do ibin=0, 400
           do idimn=1, 3
              idis_velo_solv(idimn,ibin) = 0
              i_total_v(idimn)=0
           end do
        end do
        
        n_solv=inmpc%n_all_solv
        do is=1, n_solv
           do idimn=1, 3
              ivelo=velo_solv_mpc(idimn,is)*10.0+200
              !!ivelo=(velo_solv_mpc(idimn,is)+200)*10.0
              if(ivelo >= 0 .and. ivelo <= 400) then
                 idis_velo_solv(idimn,ivelo)=idis_velo_solv(idimn,ivelo) + 1
                 i_total_v(idimn)=i_total_v(idimn)+1
              end if
           end do
        end do
        do iii=0, 400
           write(*,*) "dis x_velo ",istep, (iii-200)*1.0/10.0+1.0/20.0, idis_velo_solv(1,iii)*1.0/i_total_v(1) 
           write(*,*) "dis y_velo ",istep, (iii-200)*1.0/10.0+1.0/20.0, idis_velo_solv(2,iii)*1.0/i_total_v(2) 
           write(*,*) "dis z_velo ",istep, (iii-200)*1.0/10.0+1.0/20.0, idis_velo_solv(3,iii)*1.0/i_total_v(3) 
        end do
     end if
  end if

end subroutine simu_collision_mpc
