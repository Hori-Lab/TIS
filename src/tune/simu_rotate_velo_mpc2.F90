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
subroutine simu_rotate_velo_mpc2(velo_mp, irep)

  use const_maxsize, only : PREC, MXSOLV_MPC, MXMP, MXGRID_N_MPC
  use const_physical, only : F_PI
  use var_setp, only : irand
  use var_struct, only : nmp_real, cmass_mp
  use var_mpc, only : inmpc, ngrid_mpc, isolv2grid_mpc, imp2grid_mpc, &
       cmass_solv_mpc, velo_solv_mpc
  
  implicit none
  ! --------------------------------------------------------------------
  real(PREC), intent(inout) :: velo_mp(:, :, :)
  integer, intent(in) :: irep

  ! --------------------------------------------------------------------
  ! function
  real(PREC) :: recipe_rand01

  ! --------------------------------------------------------------------
  !local variables
  integer :: idimn, is, imp, igr,itmp_gr
  integer :: n_solv, n_grid_all
  real(PREC) :: P_total_grid(3, 0:MXGRID_N_MPC), velo_cm_grid(3, 0:MXGRID_N_MPC)
  real(PREC) :: P_total_grid_af(3, 0:MXGRID_N_MPC)
  
  real(PREC) :: total_mass_grid(0:MXGRID_N_MPC)
  real(PREC) :: vector_n_grid(3, 0:MXGRID_N_MPC)
  integer :: icount_grid(0:MXGRID_N_MPC), total_count
  !real(PREC) :: rand_theta
  real(PREC) :: rand_phi
  real(PREC) :: zzz, szzz

  real(PREC) :: delta_velo(3)
  real(PREC) :: delta_velo_new(3)
  real(PREC) :: coef_inprdct, vector_exprdct(3),coef_A
  real(PREC) :: cos_alp, sin_alp

  ! --------------------------------------------------------------------
  !! ≤Û≈æ¶¡§Œcos, sin∑◊ªª
  cos_alp = cos(F_PI*(inmpc%rotate_angle_colli)/(180.0e0_PREC))
  sin_alp = sin(F_PI*(inmpc%rotate_angle_colli)/(180.0e0_PREC))

  
  ! --------------------------------------------------------------------
  !! Initialization
  !! P_total_grid(1:3, igr):  total momentum vector in grid [igr]
  !! velo_cm_grid(1:3, igr):  velocity of center of mass in grid [igr]
  !! total_mass_grid(igr):  total mass in grid [igr]
  !! icount_grid(igr):  particle (solvent, system) number in grid [igr]
  !! vector_n_grid(1:3, igr): rotation axis vector in grid [igr]
  
  n_grid_all = ngrid_mpc(1)*ngrid_mpc(2)*ngrid_mpc(3) !! total grid number
  do igr = 0, n_grid_all
     P_total_grid(1:3, igr) = 0.0e0_PREC
     P_total_grid_af(1:3, igr) = 0.0e0_PREC

     velo_cm_grid(1:3, igr) = 0.0e0_PREC
     total_mass_grid(igr) = 0.0e0_PREC
     icount_grid(igr) = 0
     vector_n_grid(1:3, igr) = 0.0e0_PREC
  end do
  
  ! --------------------------------------------------------------------
  !! Calculation of velo_cm_grid(igr)  &  vector_n_grid(igr)
  !! P_total_grid(1:3, igr):  total momentum vector in grid [igr]
  !! velo_cm_grid(1:3, igr):  velocity of center of mass in grid [igr]
  !! vector_n_grid(1:3, igr): rotation axis vector in grid [igr]
  
  !!(pre-calculation)
  !! sum  part of solvent particle 
  n_solv = inmpc%n_all_solv
  do is = 1, n_solv
     itmp_gr = isolv2grid_mpc(is)
     P_total_grid(1:3, itmp_gr) = P_total_grid(1:3, itmp_gr) &
          + cmass_solv_mpc(is) * velo_solv_mpc(1:3, is)
     total_mass_grid(itmp_gr) = total_mass_grid(itmp_gr) + cmass_solv_mpc(is)
     icount_grid(itmp_gr) = icount_grid(itmp_gr) + 1
  end do
  
  !!(pre-calculation)
  !! sum part of system particle
  do imp = 1, nmp_real
     itmp_gr = imp2grid_mpc(imp)
     P_total_grid(1:3, itmp_gr)= P_total_grid(1:3, itmp_gr) &
          + cmass_mp(imp) * velo_mp(1:3, imp, irep)
     total_mass_grid(itmp_gr) = total_mass_grid(itmp_gr) + cmass_mp(imp)
     icount_grid(itmp_gr) = icount_grid(itmp_gr) + 1
  end do
  
  !!(calculation velocity of center of mass)
  total_count = 0
  do igr = 0, n_grid_all - 1
     if(icount_grid(igr) > 0) then
        total_count = total_count + icount_grid(igr)
        !! velocity of center of mass in grid [igr]
        velo_cm_grid(1:3, igr) = P_total_grid(1:3, igr) / total_mass_grid(igr)
        
        !! randam axis for rotation in grid [igr]
        !! rand_theta[0:pi], rand_phi[0:2pi]
        !! rand_theta = F_PI*recipe_rand01(irand)

        !! randam variable zzz[ -1 : 1]
        zzz = 2.0e0_PREC * recipe_rand01(irand) - 1.0e0_PREC
        !rand_theta = acos(zzz)
        szzz = 1.0 - zzz**2
        if(szzz < 0) then
           szzz = 0
        else
           szzz = sqrt(szzz)
        end if
        rand_phi = 2.0e0_PREC * F_PI * recipe_rand01(irand)

        !vector_n_grid(1, igr) = sin(rand_theta)*cos(rand_phi) !!n_x
        !vector_n_grid(2, igr) = sin(rand_theta)*sin(rand_phi) !!n_y
        !vector_n_grid(3, igr) = cos(rand_theta)               !!n_z

        vector_n_grid(1, igr) = szzz*cos(rand_phi)
        vector_n_grid(2, igr) = szzz*sin(rand_phi)
        vector_n_grid(3, igr) = zzz

     else
        !! nothing 
     endif
  end do
  !!write(*,*) "check particle number in total grid",total_count, nmp+n_solv


  !*****************************************************
  !*****************************************************
  !! •È•Û•¿•‡§ º¥(nx,ny,nz)§Ú√Êø¥§À ¶¡≤Û≈æ

  !!rotation <solvent particel in mpc>
  n_solv = inmpc%n_all_solv
  do is = 1, n_solv
     itmp_gr = isolv2grid_mpc(is)
     coef_inprdct = 0.0e0_PREC
     do idimn = 1, 3
        !!¶§v_0•Ÿ•Ø•»•Î∫Ó¿Æ
        delta_velo(idimn) = velo_solv_mpc(idimn, is) - velo_cm_grid(idimn, itmp_gr)
        !! ∆‚¿—∑◊ªª (n*¶§v_0)
        coef_inprdct = coef_inprdct + vector_n_grid(idimn, itmp_gr) * delta_velo(idimn)
     end do
     !! ≥∞¿—•Ÿ•Ø•»•Î§Œ∫Ó¿Æ (n°ﬂ¶§v)
     vector_exprdct(1) = vector_n_grid(2, itmp_gr)*delta_velo(3)-vector_n_grid(3, itmp_gr)*delta_velo(2)
     vector_exprdct(2) = vector_n_grid(3, itmp_gr)*delta_velo(1)-vector_n_grid(1, itmp_gr)*delta_velo(3)
     vector_exprdct(3) = vector_n_grid(1, itmp_gr)*delta_velo(2)-vector_n_grid(2, itmp_gr)*delta_velo(1)
     
     coef_A = (1.0e0_PREC - cos_alp)*coef_inprdct
     do idimn = 1, 3
        !!¶¡≤Û≈æ∏Â ¶§v_1 •Ÿ•Ø•»•Î 
        delta_velo_new(idimn) = coef_A * vector_n_grid(idimn, itmp_gr) + cos_alp * delta_velo(idimn) + sin_alp * vector_exprdct(idimn)
        
        !! ¶¡≤Û≈æ∏Â solvent velocity
        velo_solv_mpc(idimn, is) = velo_cm_grid(idimn, itmp_gr) + delta_velo_new(idimn)
     end do

  end do
  
  !!rotation <system particel in mpc>
  do imp=1, nmp_real

     itmp_gr=imp2grid_mpc(imp)
     coef_inprdct=0.0e0_PREC
     do idimn=1, 3
        !!¶§v_0•Ÿ•Ø•»•Î∫Ó¿Æ
        delta_velo(idimn) = velo_mp(idimn, imp, irep) - velo_cm_grid(idimn, itmp_gr)
        !! ∆‚¿—∑◊ªª (n*¶§v_0)
        coef_inprdct = coef_inprdct + vector_n_grid(idimn, itmp_gr) * delta_velo(idimn)
     end do
     !! ≥∞¿—•Ÿ•Ø•»•Î§Œ∫Ó¿Æ (n°ﬂ¶§v)
     vector_exprdct(1) = vector_n_grid(2, itmp_gr)*delta_velo(3)-vector_n_grid(3, itmp_gr)*delta_velo(2)
     vector_exprdct(2) = vector_n_grid(3, itmp_gr)*delta_velo(1)-vector_n_grid(1, itmp_gr)*delta_velo(3)
     vector_exprdct(3) = vector_n_grid(1, itmp_gr)*delta_velo(2)-vector_n_grid(2, itmp_gr)*delta_velo(1)

     coef_A=(1.0e0_PREC - cos_alp)*coef_inprdct
     do idimn=1, 3
        !!¶¡≤Û≈æ∏Â ¶§v_1 •Ÿ•Ø•»•Î 
        delta_velo_new(idimn) = coef_A * vector_n_grid(idimn, itmp_gr) + cos_alp * delta_velo(idimn) + sin_alp * vector_exprdct(idimn)
        
        !! ¶¡≤Û≈æ∏Â solvent velocity
        velo_mp(idimn, imp, irep) = velo_cm_grid(idimn, itmp_gr) + delta_velo_new(idimn)
     end do
  end do


end subroutine simu_rotate_velo_mpc2
