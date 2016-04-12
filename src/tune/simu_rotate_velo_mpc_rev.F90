!!**********************************************
!! input 
!!integer:: n_solv

!!integer, save :: isolv2grid_mpc(MXSOLV_MPC) !! mpc solvent particles devied into grid
!!real(PREC), save :: velo_solv_mpc(3, MXSOLV_MPC) ! velocity of solvent particle (MPC)
!!real(PREC), save :: cmass_solv_mpc(MXSOLV_MPC) ! velocity of solvent particle (MPC)

!!integer:: nmp
!!integer, save :: imp2grid_mpc(MXMP) !! system particles devied into grid
!!real(PREC), intent(inout) :: velo_mp(3, MXMP) !! velocity of system particle
!!use var_struct, only : cmass_mp

!!
!! output
!!real(PREC), save :: velo_solv_mpc(3, MXSOLV_MPC) ! velocity of solvent particle (MPC)
!!real(PREC), intent(inout) :: velo_mp(3, MXMP) !! velocity of system particle
!!**********************************************

!! rotate velocity vector (solvent, system particle) at collision step in mpc dynamics
subroutine simu_rotate_velo_mpc_rev(velo_mp, irep)

  use const_maxsize, only : PREC, MXSOLV_MPC, MXMP, MXGRID_N_MPC
  use const_physical, only : F_PI
  use var_setp, only : irand
  use var_struct, only : nmp_real, cmass_mp
  use var_mpc, only: inmpc, ngrid_mpc,isolv2grid_mpc,imp2grid_mpc,&
       cmass_solv_mpc,velo_solv_mpc
  
  implicit none
  !!-------------------------------------------
  real(PREC), intent(inout) :: velo_mp(:, :, :)
  integer, intent(in) :: irep
  !!-------------------------------------------
  ! --------------------------------------------------------------------
  ! function
  real(PREC) :: recipe_rand01
  ! --------------------------------------------------------------------

  !-----------------
  !local variables
  integer :: idimn, jdimn, is, imp, igr,itmp_gr
  integer :: n_solv, n_grid_all
  real(PREC) :: P_total_grid(3, 0:MXGRID_N_MPC), velo_cm_grid(3, 0:MXGRID_N_MPC)
  real(PREC) :: P_total_grid_af(3, 0:MXGRID_N_MPC)
  !real(PREC) :: P_total_grid_bf_gizi(3, 0:MXGRID_N_MPC)
  
  real(PREC) :: E_total_grid(0:MXGRID_N_MPC)
  real(PREC) :: E_total_grid_af(0:MXGRID_N_MPC)

  real(PREC) :: total_mass_grid(0:MXGRID_N_MPC)
  real(PREC) :: vector_n_grid(3, 0:MXGRID_N_MPC)
  integer :: icount_grid(0:MXGRID_N_MPC), total_count
  real(PREC) :: rand_theta, rand_phi
  real(PREC) :: zzz
  !real(PREC) :: delta_velo(3)
  !real(PREC) :: delta_velo_new(3)
  !real(PREC) :: coef_inprdct, vector_exprdct(3), coef_A
  real(PREC) :: cos_alp, sin_alp
  !real(PREC) :: check_velo2_bf, check_velo2_af
  real(PREC) :: check_total_p_bf(3),check_total_p_af(3),E_total_bff,E_total_aff

!!===new variabel==================-
  real(PREC) :: Rotate_matrix_grid(3, 3, 0:MXGRID_N_MPC)
  real(PREC) :: velo_cm_rev_grid(3, 0:MXGRID_N_MPC)
  real(PREC) :: cos_alp_m  
  real(PREC) :: vect_nxy_grid, vect_nxz_grid, vect_nyz_grid
  real(PREC) :: vect_nxx_grid, vect_nyy_grid, vect_nzz_grid
  real(PREC) :: tmp_velo(3)
!!===new variabel==================-
  
  !! ≤Û≈æ¶¡§Œcos, sin∑◊ªª
  cos_alp = cos(F_PI*(inmpc%rotate_angle_colli)/(180.0e0_PREC))
  sin_alp = sin(F_PI*(inmpc%rotate_angle_colli)/(180.0e0_PREC))
  cos_alp_m = (1.0e0_PREC - cos_alp)

  
!!*************************************************
  !! Initialization
  !! P_total_grid(1:3, igr):  total momentum vector in grid [igr]
  !! velo_cm_grid(1:3, igr):  velocity of center of mass in grid [igr]
  !! total_mass_grid(igr):  total mass in grid [igr]
  !! icount_grid(igr):  particle (solvent, system) number in grid [igr]
  !! vector_n_grid(1:3, igr): rotation axis vector in grid [igr]
  
  n_grid_all=ngrid_mpc(1)*ngrid_mpc(2)*ngrid_mpc(3) !! total grid number
  !!write(*,*) "pass kanada1", n_grid_all
  !!write(*,*) "MXGRID_N_MPC=",MXGRID_N_MPC
  do igr=0, n_grid_all
     P_total_grid(1:3, igr)=0.0e0_PREC

     P_total_grid_af(1:3, igr)=0.0e0_PREC
     !!P_total_grid_bf_gizi(1:3, igr)=0.0e0_PREC

    !! E_total_grid(igr)=0.0e0_PREC
    !! E_total_grid_af(igr)=0.0e0_PREC

     velo_cm_grid(1:3, igr) = 0.0e0_PREC
     total_mass_grid(igr)=0.0e0_PREC
     icount_grid(igr)=0
     vector_n_grid(1:3, igr)=0.0e0_PREC
   
     Rotate_matrix_grid(1:3, 1:3, igr)=0.0e0_PREC
     velo_cm_rev_grid(1:3, igr)=0.0e0_PREC     
	

  end do
!!*************************************************  
  !!write(*,*) "pass kanada2"
  
!!*************************************************
  !! Calculation of velo_cm_grid(igr)  &  vector_n_grid(igr)
  !! P_total_grid(1:3, igr):  total momentum vector in grid [igr]
  !! velo_cm_grid(1:3, igr):  velocity of center of mass in grid [igr]
  !! vector_n_grid(1:3, igr): rotation axis vector in grid [igr]
  
  !!(pre-calculation)
  !! sum  part of solvent particle 
  n_solv=inmpc%n_all_solv
  do is=1, n_solv
     itmp_gr=isolv2grid_mpc(is)
     do idimn=1, 3
        P_total_grid(idimn, itmp_gr)= P_total_grid(idimn, itmp_gr) + cmass_solv_mpc(is) * velo_solv_mpc(idimn, is)
        !!P_total_grid_bf_gizi(idimn, itmp_gr)= P_total_grid_bf_gizi(idimn, itmp_gr) + abs(cmass_solv_mpc(is) * velo_solv_mpc(idimn, is))
        !!if(inmpc%i_flag_check_mpc ==1)then
        !!   E_total_grid(itmp_gr)= E_total_grid(itmp_gr) + 1.0e0_PREC/2.0e0_PREC*cmass_solv_mpc(is) * velo_solv_mpc(idimn, is)* velo_solv_mpc(idimn, is)
        !!endif
     end do
     total_mass_grid(itmp_gr)=total_mass_grid(itmp_gr) + cmass_solv_mpc(is)
     icount_grid(itmp_gr)=icount_grid(itmp_gr)+1
  end do
  !!(pre-calculation)
  !! sum part of system particle
  do imp=1, nmp_real
     itmp_gr=imp2grid_mpc(imp)
     do idimn=1, 3
        P_total_grid(idimn, itmp_gr)= P_total_grid(idimn, itmp_gr) + cmass_mp(imp) * velo_mp(idimn, imp, irep)
        
        !!P_total_grid_bf_gizi(idimn, itmp_gr)= P_total_grid_bf_gizi(idimn, itmp_gr) + abs(cmass_mp(imp) * velo_mp(idimn, imp))
        !!if(inmpc%i_flag_check_mpc ==1)then
        !!   E_total_grid(itmp_gr)= E_total_grid(itmp_gr) + 1.0e0_PREC/2.0e0_PREC*cmass_mp(imp) * velo_mp(idimn, imp)* velo_mp(idimn, imp)
        !!endif
     end do
     total_mass_grid(itmp_gr)=total_mass_grid(itmp_gr) + cmass_mp(imp)
     icount_grid(itmp_gr)=icount_grid(itmp_gr)+1
  end do
  
  
  
  !!(calculation velocity of center of mass)
  total_count=0
  do igr=0, n_grid_all-1
     if(icount_grid(igr)>0)then
        total_count=total_count+icount_grid(igr)
        !! velocity of center of mass in grid [igr]
        velo_cm_grid(1:3, igr)=P_total_grid(1:3, igr)/total_mass_grid(igr)
        
        !! randam axis for rotation in grid [igr]
        !! rand_theta[0:pi], rand_phi[0:2pi]
        !! rand_theta = F_PI*recipe_rand01(irand)

        !! randam variable zzz[ -1 : 1]
        zzz = 2.0e0_PREC * recipe_rand01(irand) - 1.0e0_PREC
        rand_theta = acos(zzz)
        !! rand_theta = 2.0e0_PREC * acos(sqrt(1.0e0_PREC - recipe_rand01(irand)))        
        rand_phi = 2.0e0_PREC*F_PI*recipe_rand01(irand)

        vector_n_grid(1, igr)=sin(rand_theta)*cos(rand_phi) !!n_x
        vector_n_grid(2, igr)=sin(rand_theta)*sin(rand_phi) !!n_y
        vector_n_grid(3, igr)=cos(rand_theta)               !!n_z

        !=========================================
        ! calculate Rotational matrix:  Rotate_matrix  for each grid
        !=========================================
        vect_nxx_grid = vector_n_grid(1, igr) * vector_n_grid(1, igr)
        vect_nyy_grid = vector_n_grid(2, igr) * vector_n_grid(2, igr)
        vect_nzz_grid = vector_n_grid(3, igr) * vector_n_grid(3, igr)

        vect_nxy_grid = vector_n_grid(1, igr) * vector_n_grid(2, igr)
        vect_nxz_grid = vector_n_grid(1, igr) * vector_n_grid(3, igr)
        vect_nyz_grid = vector_n_grid(2, igr) * vector_n_grid(3, igr)

        !! diagonal part
          Rotate_matrix_grid(1, 1, igr) = cos_alp_m * vect_nxx_grid + cos_alp
          Rotate_matrix_grid(2, 2, igr) = cos_alp_m * vect_nyy_grid + cos_alp
          Rotate_matrix_grid(3, 3, igr) = cos_alp_m * vect_nzz_grid + cos_alp

       !! off-diagonal part
          Rotate_matrix_grid(1, 2, igr) = cos_alp_m * vect_nxy_grid - sin_alp * vector_n_grid(3, igr)
          Rotate_matrix_grid(1, 3, igr) = cos_alp_m * vect_nxz_grid + sin_alp * vector_n_grid(2, igr)

          Rotate_matrix_grid(2, 1, igr) = cos_alp_m * vect_nxy_grid + sin_alp * vector_n_grid(3, igr)
          Rotate_matrix_grid(2, 3, igr) = cos_alp_m * vect_nyz_grid - sin_alp * vector_n_grid(1, igr)

          Rotate_matrix_grid(3, 1, igr) = cos_alp_m * vect_nxz_grid - sin_alp * vector_n_grid(2, igr)
          Rotate_matrix_grid(3, 2, igr) = cos_alp_m * vect_nyz_grid + sin_alp * vector_n_grid(1, igr)


        !===========================================
        ! calculate  vector [V_cm-R*V_cm]: velo_cm_rev_grid(3, 0:MXGRID_N_MPC)
        !===========================================
        do idimn=1, 3
            velo_cm_rev_grid(idimn, igr) = velo_cm_grid(idimn, igr)
           do jdimn=1, 3
            velo_cm_rev_grid(idimn, igr) = velo_cm_rev_grid(idimn, igr) - Rotate_matrix_grid(idimn, jdimn, igr) * velo_cm_grid(jdimn, igr)
           end do
        end do


     else
        !! nothing 
     endif
  end do
  !!write(*,*) "check particle number in total grid",total_count, nmp+n_solv
  
!*****************************************************
!*****************************************************
!! •È•Û•¿•‡§ º¥(nx,ny,nz)§Ú√Êø¥§À ¶¡≤Û≈æ

  !!rotation <solvent particel in mpc>
  n_solv=inmpc%n_all_solv
  do is=1, n_solv
     itmp_gr=isolv2grid_mpc(is)
     do idimn=1, 3
        tmp_velo(idimn) = velo_cm_rev_grid(idimn, itmp_gr)
        do jdimn=1, 3
           tmp_velo(idimn) = tmp_velo(idimn) + Rotate_matrix_grid(idimn, jdimn, itmp_gr) * velo_solv_mpc(jdimn, is)
        end do
     end do
     velo_solv_mpc(1:3, is) = tmp_velo(1:3)


     !!coef_inprdct=0.0e0_PREC
     !!do idimn=1, 3
     !!   !!¶§v_0•Ÿ•Ø•»•Î∫Ó¿Æ
     !!   delta_velo(idimn) = velo_solv_mpc(idimn, is) - velo_cm_grid(idimn, itmp_gr)
     !!   !! ∆‚¿—∑◊ªª (n*¶§v_0)
     !!   coef_inprdct = coef_inprdct + vector_n_grid(idimn, itmp_gr) * delta_velo(idimn)
     !!end do
     !! ≥∞¿—•Ÿ•Ø•»•Î§Œ∫Ó¿Æ (n°ﬂ¶§v)
     !!vector_exprdct(1) = vector_n_grid(2, itmp_gr)*delta_velo(3)-vector_n_grid(3, itmp_gr)*delta_velo(2)
     !!vector_exprdct(2) = vector_n_grid(3, itmp_gr)*delta_velo(1)-vector_n_grid(1, itmp_gr)*delta_velo(3)
     !!vector_exprdct(3) = vector_n_grid(1, itmp_gr)*delta_velo(2)-vector_n_grid(2, itmp_gr)*delta_velo(1)
     !!
     !!coef_A=(1.0e0_PREC - cos_alp)*coef_inprdct
     !!do idimn=1, 3
     !!   !!¶¡≤Û≈æ∏Â ¶§v_1 •Ÿ•Ø•»•Î 
     !!   delta_velo_new(idimn) = coef_A * vector_n_grid(idimn, itmp_gr) + cos_alp * delta_velo(idimn) + sin_alp * vector_exprdct(idimn)
     !!   
     !!   !! ¶¡≤Û≈æ∏Â solvent velocity
     !!   velo_solv_mpc(idimn, is) = velo_cm_grid(idimn, itmp_gr) + delta_velo_new(idimn)
     !!end do
     !!
     !!check_velo2_bf=0.0e0_PREC
     !!check_velo2_af=0.0e0_PREC
     !!do idimn=1, 3
     !!   check_velo2_bf=check_velo2_bf+delta_velo(idimn)*delta_velo(idimn)
     !!   check_velo2_af=check_velo2_af+delta_velo_new(idimn)*delta_velo_new(idimn)
     !!end do
     !!write(*,*) "cehck velo2_bf,velo2_af",check_velo2_bf- check_velo2_af
  end do
  
  !!rotation <system particel in mpc>
  do imp=1, nmp_real
     
     !!check_velo2_bf=0.0e0_PREC
     !!do idimn=1, 3
     !!   check_velo2_bf=check_velo2_bf+velo_mp(idimn,imp)*velo_mp(idimn,imp)
     !!end do

     itmp_gr=imp2grid_mpc(imp)
     do idimn=1, 3
        tmp_velo(idimn) = velo_cm_rev_grid(idimn, itmp_gr)
         do jdimn=1, 3
            tmp_velo(idimn) = tmp_velo(idimn) + Rotate_matrix_grid(idimn, jdimn, itmp_gr) *  velo_mp(jdimn, imp, irep)
         end do
     end do
     velo_mp(1:3, imp, irep) = tmp_velo(1:3)



     !!coef_inprdct=0.0e0_PREC
     !!do idimn=1, 3
     !!   !!¶§v_0•Ÿ•Ø•»•Î∫Ó¿Æ
     !!   delta_velo(idimn) = velo_mp(idimn, imp) - velo_cm_grid(idimn, itmp_gr)
     !!   !! ∆‚¿—∑◊ªª (n*¶§v_0)
     !!   coef_inprdct = coef_inprdct + vector_n_grid(idimn, itmp_gr) * delta_velo(idimn)
     !!end do
     !! ≥∞¿—•Ÿ•Ø•»•Î§Œ∫Ó¿Æ (n°ﬂ¶§v)
     !!vector_exprdct(1) = vector_n_grid(2, itmp_gr)*delta_velo(3)-vector_n_grid(3, itmp_gr)*delta_velo(2)
     !!vector_exprdct(2) = vector_n_grid(3, itmp_gr)*delta_velo(1)-vector_n_grid(1, itmp_gr)*delta_velo(3)
     !!vector_exprdct(3) = vector_n_grid(1, itmp_gr)*delta_velo(2)-vector_n_grid(2, itmp_gr)*delta_velo(1)
     !!
     !!coef_A=(1.0e0_PREC - cos_alp)*coef_inprdct
     !!do idimn=1, 3
     !!   !!¶¡≤Û≈æ∏Â ¶§v_1 •Ÿ•Ø•»•Î 
     !!   delta_velo_new(idimn) = coef_A * vector_n_grid(idimn, itmp_gr) + cos_alp * delta_velo(idimn) + sin_alp * vector_exprdct(idimn)
     !!   
     !!   !! ¶¡≤Û≈æ∏Â solvent velocity
     !!   velo_mp(idimn, imp) = velo_cm_grid(idimn, itmp_gr) + delta_velo_new(idimn)
     !!end do
     
     !!check_velo2_bf=0.0e0_PREC
     !!check_velo2_af=0.0e0_PREC
     !!do idimn=1, 3
     !!   check_velo2_bf=check_velo2_bf+delta_velo(idimn)*delta_velo(idimn)
     !!   check_velo2_af=check_velo2_af+delta_velo_new(idimn)*delta_velo_new(idimn)
     !!end do
     !!write(*,*) "cehck velo2_bf,velo2_af",check_velo2_bf, check_velo2_af,check_velo2_bf-check_velo2_af
     
     !!check_velo2_bf=0.0e0_PREC
     !!check_velo2_af=0.0e0_PREC
     !!do idimn=1, 3
     !!   !!check_velo2_bf=check_velo2_bf+
     !!   check_velo2_af=check_velo2_af+velo_mp(idimn,imp)*velo_mp(idimn,imp)
     !!end do
     !!write(*,*) "cehck velo2_bf,velo2_af",check_velo2_bf, check_velo2_af,check_velo2_bf-check_velo2_af
    
  end do


  if(inmpc%i_flag_check_mpc ==1)then
     n_solv=inmpc%n_all_solv
     do is=1, n_solv
        itmp_gr=isolv2grid_mpc(is)
        do idimn=1, 3
           !!P_total_grid_af(idimn, itmp_gr)= P_total_grid_af(idimn, itmp_gr) + abs(cmass_solv_mpc(is) * velo_solv_mpc(idimn, is))
           P_total_grid_af(idimn, itmp_gr)= P_total_grid_af(idimn, itmp_gr) + (cmass_solv_mpc(is) * velo_solv_mpc(idimn, is))
           E_total_grid_af(itmp_gr)= E_total_grid_af(itmp_gr) + 1.0e0_PREC/2.0e0_PREC*cmass_solv_mpc(is) * velo_solv_mpc(idimn, is)* velo_solv_mpc(idimn, is)
        end do
     end do
     do imp=1, nmp_real
        itmp_gr=imp2grid_mpc(imp)
        do idimn=1, 3
           !!P_total_grid_af(idimn, itmp_gr)= P_total_grid_af(idimn, itmp_gr) + abs(cmass_mp(imp) * velo_mp(idimn, imp))
           P_total_grid_af(idimn, itmp_gr)= P_total_grid_af(idimn, itmp_gr) + (cmass_mp(imp) * velo_mp(idimn, imp, irep))
           E_total_grid_af(itmp_gr)= E_total_grid_af(itmp_gr) + 1.0e0_PREC/2.0e0_PREC*cmass_mp(imp) * velo_mp(idimn, imp, irep)* velo_mp(idimn, imp, irep)
        end do
     end do
     
     !!  do igr=0, n_grid_all
     !!        write(*,*) "New check P_bf",P_total_grid(1, igr),P_total_grid(2, igr),P_total_grid(3, igr)
     !!        !!   write(*,*) "New check P_gf",P_total_grid_bf_gizi(1, igr),P_total_grid_bf_gizi(2, igr),P_total_grid_bf_gizi(3, igr)
     !!        write(*,*) "New check P_af",P_total_grid_af(1, igr),P_total_grid_af(2, igr),P_total_grid_af(3, igr)
     !!        !!   write(*,*) "New check E_bf,af",E_total_grid(igr),E_total_grid_af(igr)
     !!  end do
     

     
     check_total_p_bf(1)=0.0e0_PREC
     check_total_p_bf(2)=0.0e0_PREC
     check_total_p_bf(3)=0.0e0_PREC
     check_total_p_af(1)=0.0e0_PREC
     check_total_p_af(2)=0.0e0_PREC
     check_total_p_af(3)=0.0e0_PREC
     E_total_bff=0.0e0_PREC
     E_total_aff=0.0e0_PREC
  
     do igr=0, n_grid_all
        check_total_p_bf(1)=check_total_p_bf(1)+P_total_grid(1, igr)
        check_total_p_bf(2)=check_total_p_bf(2)+P_total_grid(2, igr)
        check_total_p_bf(3)=check_total_p_bf(3)+P_total_grid(3, igr)
        !     check_total_p_bf(1)=check_total_p_bf(1)+P_total_grid_bf_gizi(1, igr)
        !     check_total_p_bf(2)=check_total_p_bf(2)+P_total_grid_bf_gizi(2, igr)
        !     check_total_p_bf(3)=check_total_p_bf(3)+P_total_grid_bf_gizi(3, igr)
        check_total_p_af(1)=check_total_p_af(1)+P_total_grid_af(1, igr)
        check_total_p_af(2)=check_total_p_af(2)+P_total_grid_af(2, igr)
        check_total_p_af(3)=check_total_p_af(3)+P_total_grid_af(3, igr)
        E_total_bff=E_total_bff+E_total_grid(igr)
        E_total_aff=E_total_aff+E_total_grid_af(igr)
     end do
     write(*,*) "Tot check P_x_bf,af",check_total_p_bf(1),check_total_p_af(1)
     write(*,*) "Tot check P_y_bf,af",check_total_p_bf(2),check_total_p_af(2)
     write(*,*) "Tot check P_z_bf,af",check_total_p_bf(3),check_total_p_af(3)
     write(*,*) "Tot check E_bf,af",E_total_bff,E_total_aff
  end if


!!  check_total_p_af(1)=0.0e0_PREC
!!  check_total_p_af(2)=0.0e0_PREC
!!  check_total_p_af(3)=0.0e0_PREC
  !#*********************

!!  n_solv=inmpc%n_all_solv
!!  do is=1, n_solv
!!     check_total_p_af(1) = check_total_p_af(1) + cmass_solv_mpc(is) * velo_solv_mpc(1, is)
!!     check_total_p_af(2) = check_total_p_af(2) + cmass_solv_mpc(is) * velo_solv_mpc(2, is)
!!     check_total_p_af(3) = check_total_p_af(3) + cmass_solv_mpc(is) * velo_solv_mpc(3, is)
!!  end do

  !!(pre-calculation)
  !! sum part of system particle
!!  do imp=1, nmp
!!     check_total_p_af(1) = check_total_p_af(1) + cmass_mp(imp) * velo_mp(1, imp)
!!     check_total_p_af(2) = check_total_p_af(2) + cmass_mp(imp) * velo_mp(2, imp)
!!     check_total_p_af(3) = check_total_p_af(3) + cmass_mp(imp) * velo_mp(3, imp)
!!  end do
  !!#********************
!!  write(*,*) "Tot check Paf2_x,y,zf",check_total_p_af(1),check_total_p_af(2),check_total_p_af(3)


end subroutine simu_rotate_velo_mpc_rev
