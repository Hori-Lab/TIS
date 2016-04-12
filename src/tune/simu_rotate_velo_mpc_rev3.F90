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
subroutine simu_rotate_velo_mpc_rev3(velo_mp, irep, velo_solv_l, velo_mp_l)

  use const_maxsize, only : PREC, MXSOLV_MPC, MXMP, MXGRID_N_MPC
  use const_physical, only : F_PI
  use var_setp, only : irand
  use var_struct, only : nmp_real, cmass_mp
  use var_mpc, only: inmpc, ngrid_mpc, cmass_solv_mpc, velo_solv_mpc, &
                     igrid_process, lgrid_l, igrid_l, &
                     nis_l, nimp_l, isolv2grid_mpc_l, imp2grid_mpc_l
  use mpiconst
  
  implicit none
  ! --------------------------------------------------------------------
  real(PREC), intent(inout) :: velo_mp(:, :, :)
  real(PREC), intent(inout) :: velo_solv_l(:, :), velo_mp_l(:, :, :)
  integer, intent(in) :: irep

  ! --------------------------------------------------------------------
  ! function
  real(PREC) :: recipe_rand01
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !local variables
  integer :: idimn, jdimn, is, imp, igr, itmp_gr
  integer :: is_l, imp_l
  integer :: n_grid_all
  real(PREC) :: P_total_grid(3, 0:MXGRID_N_MPC)
  real(PREC) :: velo_cm_grid(3, 0:MXGRID_N_MPC)
  real(PREC) :: total_mass_grid(0:MXGRID_N_MPC)
  real(PREC) :: vector_n_grid(3, 0:MXGRID_N_MPC)
  integer :: icount_grid(0:MXGRID_N_MPC), total_count
  !real(PREC) :: rand_theta
  real(PREC) :: rand_phi
  real(PREC) :: zzz, szzz, dzzz
  real(PREC) :: cos_alp, sin_alp

  !!===new variabel==================-
  real(PREC) :: Rotate_matrix_grid(3, 3, 0:MXGRID_N_MPC)
  real(PREC) :: velo_cm_rev_grid(3, 0:MXGRID_N_MPC)
  real(PREC) :: cos_alp_m  
  real(PREC) :: vect_nxy_grid, vect_nxz_grid, vect_nyz_grid
  real(PREC) :: vect_nxx_grid, vect_nyy_grid, vect_nzz_grid
  real(PREC) :: tmp_velo(3)
  !!===new variabel==================-
  

  ! --------------------------------------------------------------------
  !! ²óÅ¾¦Á¤Îcos, sin·×»»
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
  
  n_grid_all = ngrid_mpc(1)*ngrid_mpc(2)*ngrid_mpc(3) !! total grid number
  
  P_total_grid(1:3, 0:n_grid_all) = 0.0e0_PREC
  total_mass_grid(0:n_grid_all) = 0.0e0_PREC
  icount_grid(0:n_grid_all) = 0

  velo_cm_grid(1:3, 0:n_grid_all) = 0.0e0_PREC
  velo_cm_rev_grid(1:3, 0:n_grid_all) = 0.0e0_PREC
  vector_n_grid(1:3, 0:n_grid_all) = 0.0e0_PREC
  Rotate_matrix_grid(1:3, 1:3, 0:n_grid_all) = 0.0e0_PREC
  
  !!(pre-calculation)
  !! sum  part of solvent particle
  do is_l = 1, nis_l
     itmp_gr = isolv2grid_mpc_l(1, is_l)
     is = isolv2grid_mpc_l(2, is_l)
     P_total_grid(1:3, itmp_gr) = P_total_grid(1:3, itmp_gr) &
          + cmass_solv_mpc(is) * velo_solv_mpc(1:3, is)
     total_mass_grid(itmp_gr) = total_mass_grid(itmp_gr) + cmass_solv_mpc(is)
     icount_grid(itmp_gr) = icount_grid(itmp_gr) + 1
  end do
  
  !!(pre-calculation)
  !! sum part of system particle
  do imp_l = 1, nimp_l
     itmp_gr = imp2grid_mpc_l(1, imp_l)
     imp = imp2grid_mpc_l(2, imp_l)
     P_total_grid(1:3, itmp_gr) = P_total_grid(1:3, itmp_gr) + cmass_mp(imp) * velo_mp(1:3, imp, irep)
     total_mass_grid(itmp_gr) = total_mass_grid(itmp_gr) + cmass_mp(imp)
     icount_grid(itmp_gr) = icount_grid(itmp_gr) + 1
  end do


  ! --------------------------------------------------------------------
  !!*************************************************
  !! Calculation of velo_cm_grid(igr)  &  vector_n_grid(igr)
  !! P_total_grid(1:3, igr):  total momentum vector in grid [igr]
  !! velo_cm_grid(1:3, igr):  velocity of center of mass in grid [igr]
  !! vector_n_grid(1:3, igr): rotation axis vector in grid [igr]
  
  ! --------------------------------------------------------------------
  !!(calculation velocity of center of mass)
  total_count = 0
  dzzz = 0
!  do igr = 0, lgrid_l(myrank) - 1
  do igr = 0, n_grid_all - 1

     if(icount_grid(igr) == 0) then
	dzzz = dzzz + recipe_rand01(irand)
	dzzz = dzzz + recipe_rand01(irand)
        cycle
     end if

     total_count = total_count + icount_grid(igr)
     !! velocity of center of mass in grid [igr]
     velo_cm_grid(1:3, igr) = P_total_grid(1:3, igr) / total_mass_grid(igr)
     
     !! randam axis for rotation in grid [igr]
     !! rand_theta[0:pi], rand_phi[0:2pi]
     !! rand_theta = F_PI*recipe_rand01(irand)
     
     !! randam variable zzz[ -1 : 1]
     zzz = 2.0e0_PREC * recipe_rand01(irand) - 1.0e0_PREC
     !        rand_theta = acos(zzz)
     szzz = 1.0 - zzz**2
     if(szzz < 0) then
        szzz = 0
     else
        szzz = sqrt(szzz)
     end if
     rand_phi = 2.0e0_PREC*F_PI*recipe_rand01(irand)
     
     !        vector_n_grid(1, igr) = sin(rand_theta)*cos(rand_phi) !!n_x
     !        vector_n_grid(2, igr) = sin(rand_theta)*sin(rand_phi) !!n_y
     !        vector_n_grid(3, igr) = cos(rand_theta)               !!n_z
     
     vector_n_grid(1, igr) = szzz*cos(rand_phi)
     vector_n_grid(2, igr) = szzz*sin(rand_phi)
     vector_n_grid(3, igr) = zzz
     
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
     Rotate_matrix_grid(1, 2, igr) = cos_alp_m * vect_nxy_grid &
          - sin_alp * vector_n_grid(3, igr)
     Rotate_matrix_grid(1, 3, igr) = cos_alp_m * vect_nxz_grid &
          + sin_alp * vector_n_grid(2, igr)
     
     Rotate_matrix_grid(2, 1, igr) = cos_alp_m * vect_nxy_grid &
          + sin_alp * vector_n_grid(3, igr)
     Rotate_matrix_grid(2, 3, igr) = cos_alp_m * vect_nyz_grid &
          - sin_alp * vector_n_grid(1, igr)
     
     Rotate_matrix_grid(3, 1, igr) = cos_alp_m * vect_nxz_grid &
          - sin_alp * vector_n_grid(2, igr)
     Rotate_matrix_grid(3, 2, igr) = cos_alp_m * vect_nyz_grid &
          + sin_alp * vector_n_grid(1, igr)
     
     
     !===========================================
     ! calculate  vector [V_cm-R*V_cm]: velo_cm_rev_grid(3, 0:MXGRID_N_MPC)
     !===========================================
     do idimn = 1, 3
        velo_cm_rev_grid(idimn, igr) = velo_cm_grid(idimn, igr)
        do jdimn = 1, 3
           velo_cm_rev_grid(idimn, igr) = velo_cm_rev_grid(idimn, igr) &
                - Rotate_matrix_grid(idimn, jdimn, igr) * velo_cm_grid(jdimn, igr)
        end do
     end do
  end do

  ! --------------------------------------------------------------------
  !*****************************************************
  !*****************************************************
  !! ¥é¥ó¥À¥à¤Ê¼´(nx, ny, nz)¤òÃæ¿´¤Ë ¦Á²óÅ¾
  ! --------------------------------------------------------------------
  !!rotation <solvent particel in mpc>
  do is_l = 1, nis_l
     itmp_gr = isolv2grid_mpc_l(1, is_l)
     is = isolv2grid_mpc_l(2, is_l)
     do idimn = 1, 3
        tmp_velo(idimn) = velo_cm_rev_grid(idimn, itmp_gr)
        do jdimn = 1, 3
           tmp_velo(idimn) = tmp_velo(idimn) &
                + Rotate_matrix_grid(idimn, jdimn, itmp_gr) * velo_solv_mpc(jdimn, is)
        end do
     end do
     velo_solv_l(1:3, is) = tmp_velo(1:3)
  end do
  
  !!rotation <system particel in mpc>
  do imp_l = 1, nimp_l
     itmp_gr = imp2grid_mpc_l(1, imp_l)
     imp = imp2grid_mpc_l(2, imp_l)
     do idimn = 1, 3
        tmp_velo(idimn) = velo_cm_rev_grid(idimn, itmp_gr)
        do jdimn = 1, 3
           tmp_velo(idimn) = tmp_velo(idimn) &
                + Rotate_matrix_grid(idimn, jdimn, itmp_gr) *  velo_mp(jdimn, imp, irep)
        end do
     end do
     velo_mp_l(1:3, imp, irep) = tmp_velo(1:3)

  end do
  
end subroutine simu_rotate_velo_mpc_rev3
