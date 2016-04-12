#undef TIME
#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

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
subroutine simu_rotate_velo_mpc_rev4(velo_mp, irep)

  use const_maxsize, only : PREC, MXGRID_N_MPC
  use const_physical, only : F_PI
  use var_setp, only : mts
  use var_struct, only : nmp_real, cmass_mp
  use var_mpc, only: inmpc, ngrid_mpc, isolv2grid_mpc, imp2grid_mpc,&
                     cmass_solv_mpc, velo_solv_mpc
  use mt_stream
#ifdef TIME
  use time, only : time_s, time_e, &
                   tm_rotate_mpc01, tm_rotate_mpc02, tm_rotate_mpc03, tm_rotate_mpc04, &
                   tm_rotate_mpc05, tm_rotate_mpc06, tm_rotate_mpc07, tm_rotate_mpc08
#endif
  use mpiconst
!$ use omp_lib
  
  implicit none
  ! --------------------------------------------------------------------
  real(PREC), intent(inout) :: velo_mp(:, :, :)
  integer, intent(in) :: irep

  ! --------------------------------------------------------------------
  ! function
  !real(PREC) :: recipe_rand01
  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  !local variables
  integer :: idimn, jdimn, is, imp, igr, itmp_gr, istream
  integer :: n_solv, n_grid_all
  real(PREC) :: P_total_grid(3, 0:MXGRID_N_MPC, 0:nthreads-1)
  real(PREC) :: total_mass_grid(0:MXGRID_N_MPC, 0:nthreads-1)
  integer    :: icount_grid(0:MXGRID_N_MPC, 0:nthreads-1)
  real(PREC) :: velo_cm_grid(3, 0:MXGRID_N_MPC)
  real(PREC) :: vector_n_grid(3, 0:MXGRID_N_MPC)
  integer    :: total_count
  !real(PREC) :: rand_theta
  real(PREC) :: rand_phi
  real(PREC) :: zzz, szzz
  real(PREC) :: cos_alp, sin_alp

!!===new variabel==================-
  real(PREC) :: Rotate_matrix_grid(3, 3, 0:MXGRID_N_MPC)
  real(PREC) :: velo_cm_rev_grid(3, 0:MXGRID_N_MPC)
  real(PREC) :: cos_alp_m  
  real(PREC) :: vect_nxy_grid, vect_nxz_grid, vect_nyz_grid
  real(PREC) :: vect_nxx_grid, vect_nyy_grid, vect_nzz_grid
  real(PREC) :: tmp_velo(3)
  real(PREC) :: tmp_velo1, tmp_velo2, tmp_velo3
  real(PREC) :: velo_solv_mpc1, velo_solv_mpc2, velo_solv_mpc3
  real(PREC) :: vgrnd(2,0:MXGRID_N_MPC)
!!
  integer :: tn, nth
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

  TIME_S( tm_rotate_mpc01 )

!$omp parallel
!$omp do 
  do igr = 0, n_grid_all
     P_total_grid(1:3, igr, :) = 0.0e0_PREC

!    velo_cm_grid(1:3, igr)    = 0.0e0_PREC
     total_mass_grid(igr, :)   = 0.0e0_PREC
     icount_grid(igr, :)       = 0
!    vector_n_grid(1:3, igr)   = 0.0e0_PREC
   
!    Rotate_matrix_grid(1:3, 1:3, igr) = 0.0e0_PREC
!    velo_cm_rev_grid(1:3, igr)        = 0.0e0_PREC
  end do
!$omp end do
!$omp end parallel

  TIME_E( tm_rotate_mpc01 )

  
!!*************************************************
  !! Calculation of velo_cm_grid(igr)  &  vector_n_grid(igr)
  !! P_total_grid(1:3, igr):  total momentum vector in grid [igr]
  !! velo_cm_grid(1:3, igr):  velocity of center of mass in grid [igr]
  !! vector_n_grid(1:3, igr): rotation axis vector in grid [igr]
  
  !!(pre-calculation)
  !! sum  part of solvent particle 
  n_solv = inmpc%n_all_solv

  TIME_S( tm_rotate_mpc02 )

!$omp parallel private(tn, itmp_gr)
    tn = 0
!$  tn = omp_get_thread_num()
!$omp do
  do is = 1, n_solv
     itmp_gr = isolv2grid_mpc(is)
     P_total_grid(1:3, itmp_gr, tn) = P_total_grid(1:3, itmp_gr, tn) &
          + cmass_solv_mpc(is) * velo_solv_mpc(1:3, is)
  end do
!$omp end do
!
!$omp do
  do is = 1, n_solv
     itmp_gr = isolv2grid_mpc(is)
     total_mass_grid(itmp_gr, tn) = total_mass_grid(itmp_gr, tn) + cmass_solv_mpc(is)
     icount_grid(itmp_gr, tn) = icount_grid(itmp_gr, tn) + 1
  end do
!$omp end do
!$omp end parallel

  TIME_E( tm_rotate_mpc02 )

  !!(pre-calculation)
  !! sum part of system particle

  TIME_S( tm_rotate_mpc03 )

!!$omp parallel private(tn)
     tn = 0
!!$  tn = omp_get_thread_num()
!!$omp do
  do imp=1, nmp_real
     itmp_gr=imp2grid_mpc(imp)
     P_total_grid(1:3, itmp_gr, tn)= P_total_grid(1:3, itmp_gr, tn) + cmass_mp(imp) * velo_mp(1:3, imp, irep)
     total_mass_grid(itmp_gr, tn) = total_mass_grid(itmp_gr, tn) + cmass_mp(imp)
     icount_grid(itmp_gr, tn) = icount_grid(itmp_gr, tn) + 1
  end do
!!$omp end do
!!$omp end parallel

  TIME_E( tm_rotate_mpc03 )
  
  TIME_S( tm_rotate_mpc04 )

!$omp parallel do
  do is = 0, MXGRID_N_MPC
    do nth = 1, nthreads-1
       P_total_grid(1:3, is, 0) = P_total_grid(1:3, is, 0) &
                                + P_total_grid(1:3, is, nth)
       total_mass_grid(is, 0)   = total_mass_grid(is, 0) &
                                + total_mass_grid(is, nth)
       icount_grid(is, 0)       = icount_grid(is, 0) &
                                + icount_grid(is, nth) 
    end do
  end do
  
  TIME_E( tm_rotate_mpc04 )
  
  TIME_S( tm_rotate_mpc05 )

  vgrnd(:,:) = 0.0e0_PREC
  istream = irep
  do igr = 0, n_grid_all - 1
     if(icount_grid(igr, 0) > 0) then
!        vgrnd(1,igr) = grnd()
!        vgrnd(2,igr) = grnd()
        vgrnd(1,igr) = genrand_double1(mts(istream, 0))
        vgrnd(2,igr) = genrand_double1(mts(istream, 0))
     end if
  end do
  
  TIME_E( tm_rotate_mpc05 )

  TIME_S( tm_rotate_mpc06 )

  !!(calculation velocity of center of mass)
  total_count = 0
!$omp parallel do private(zzz,szzz,rand_phi,  &
!$omp&                    vect_nxx_grid,vect_nyy_grid,vect_nzz_grid, &
!$omp&                    vect_nxy_grid,vect_nxz_grid,vect_nyz_grid) &
!$omp& reduction(+:total_count)
  do igr = 0, n_grid_all - 1
     if(icount_grid(igr, 0) > 0) then
        total_count = total_count + icount_grid(igr, 0)
        !! velocity of center of mass in grid [igr]
        velo_cm_grid(1:3, igr) = P_total_grid(1:3, igr, 0) / total_mass_grid(igr, 0)
        
      !! randam axis for rotation in grid [igr]
      !! rand_theta[0:pi], rand_phi[0:2pi]
      !! rand_theta = F_PI*recipe_rand01(irand)

      !! randam variable zzz[ -1 : 1]
!       zzz = 2.0e0_PREC*grnd() - 1.0e0_PREC
        zzz = 2.0e0_PREC*vgrnd(1,igr) - 1.0e0_PREC
!        rand_theta = acos(zzz)
        szzz = 1.0 - zzz**2
        if(szzz < 0) then
           szzz = 0
        else
           szzz = sqrt(szzz)
        end if
!       rand_phi = 2.0e0_PREC*F_PI*grnd()
        rand_phi = 2.0e0_PREC*F_PI*vgrnd(2,igr)

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

     else
        !! nothing 
     end if
  end do
  
  TIME_E( tm_rotate_mpc06 )

!*****************************************************
!*****************************************************
!! ¥é¥ó¥À¥à¤Ê¼´(nx,ny,nz)¤òÃæ¿´¤Ë ¦Á²óÅ¾

  TIME_S( tm_rotate_mpc07 )

  !!rotation <solvent particel in mpc>
  n_solv = inmpc%n_all_solv
!$omp parallel do private(itmp_gr, tmp_velo1, tmp_velo2, tmp_velo3, &
!$omp&                    velo_solv_mpc1, velo_solv_mpc2, velo_solv_mpc3)
  do is = 1, n_solv
     itmp_gr = isolv2grid_mpc(is)

     tmp_velo1 = velo_cm_rev_grid(1, itmp_gr)
     tmp_velo2 = velo_cm_rev_grid(2, itmp_gr)
     tmp_velo3 = velo_cm_rev_grid(3, itmp_gr)

     velo_solv_mpc1 = velo_solv_mpc(1, is)
     velo_solv_mpc2 = velo_solv_mpc(2, is)
     velo_solv_mpc3 = velo_solv_mpc(3, is)

     tmp_velo1 = tmp_velo1 &
               + Rotate_matrix_grid(1, 1, itmp_gr) * velo_solv_mpc1
     tmp_velo1 = tmp_velo1 &
               + Rotate_matrix_grid(1, 2, itmp_gr) * velo_solv_mpc2
     tmp_velo1 = tmp_velo1 &
               + Rotate_matrix_grid(1, 3, itmp_gr) * velo_solv_mpc3

     tmp_velo2 = tmp_velo2 &
               + Rotate_matrix_grid(2, 1, itmp_gr) * velo_solv_mpc1
     tmp_velo2 = tmp_velo2 &
               + Rotate_matrix_grid(2, 2, itmp_gr) * velo_solv_mpc2
     tmp_velo2 = tmp_velo2 &
               + Rotate_matrix_grid(2, 3, itmp_gr) * velo_solv_mpc3

     tmp_velo3 = tmp_velo3 &
               + Rotate_matrix_grid(3, 1, itmp_gr) * velo_solv_mpc1
     tmp_velo3 = tmp_velo3 &
               + Rotate_matrix_grid(3, 2, itmp_gr) * velo_solv_mpc2
     tmp_velo3 = tmp_velo3 &
               + Rotate_matrix_grid(3, 3, itmp_gr) * velo_solv_mpc3

     velo_solv_mpc(1, is) = tmp_velo1
     velo_solv_mpc(2, is) = tmp_velo2
     velo_solv_mpc(3, is) = tmp_velo3

  end do

  TIME_E( tm_rotate_mpc07 )
  
  TIME_S( tm_rotate_mpc08 )

  !!rotation <system particel in mpc>
!$omp parallel do private(itmp_gr, tmp_velo)
  do imp = 1, nmp_real

     itmp_gr = imp2grid_mpc(imp)
     do idimn = 1, 3
        tmp_velo(idimn) = velo_cm_rev_grid(idimn, itmp_gr)
        do jdimn = 1, 3
           tmp_velo(idimn) = tmp_velo(idimn) &
                + Rotate_matrix_grid(idimn, jdimn, itmp_gr) *  velo_mp(jdimn, imp, irep)
        end do
     end do
     velo_mp(1:3, imp, irep) = tmp_velo(1:3)

  end do
  
  TIME_E( tm_rotate_mpc08 )

end subroutine simu_rotate_velo_mpc_rev4
