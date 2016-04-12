!simu_collision_mpc3
!> @brief

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

subroutine simu_collision_mpc3(irep)

  use const_maxsize
  use const_physical
!  use const_index
  use if_mloop
!  use if_write
!  use if_energy
  use var_struct, only : nmp_real, xyz_mp_rep, cmass_mp
  use var_simu, only : istep, tstep, velo_mp, tempk
  use var_mpc, only : inmpc, xyz_solv_mpc, velo_solv_mpc, ngrid_mpc, &
                      igrid_process, lgrid_l, igrid_l
  use var_replica, only : n_replica_all, n_replica_mpi

  use time, only : tm_grid_mpc, tm_rotate_mpc, tm_velo_mpc, &
                   time_s, time_e
  use mpiconst

  implicit none

  ! -----------------------------------------------------------------
  integer, intent(in) :: irep

  ! -----------------------------------------------------------------
  ! local variables
  integer :: is, n_solv, ig, ig_l, ig_p
  integer :: np, mr, ngx, ngy, ngz
  integer :: i_first = 0
  real(PREC) :: tstep_colli
  real(PREC), save :: velo_solv_l(3, MXSOLV_MPC)
  real(PREC), allocatable, save :: velo_mp_l(:, :, :)

  ! -----------------------------------------------------------------
  if(mod(istep, inmpc%nratio_colli_step) /= 0) then
     return
  end if
  
  if(i_first == 0) then
     i_first = 1

     allocate( velo_mp_l(SPACE_DIM, nmp_real, n_replica_mpi) )

     allocate( igrid_process(0:inmpc%ngrid_x, 0:inmpc%ngrid_y, 0:inmpc%ngrid_z) )
     allocate( lgrid_l(0:npar_mpi-1) )
     allocate( igrid_l(0:inmpc%ngrid_x, 0:inmpc%ngrid_y, 0:inmpc%ngrid_z) )
     lgrid_l(0:npar_mpi-1) = 0
     do ngx = 0, inmpc%ngrid_x - 1
        do ngy = 0, inmpc%ngrid_y - 1
           do ngz = 0, inmpc%ngrid_z - 1
              ig = ngrid_mpc(1)*ngrid_mpc(2)*ngz &
                   + ngrid_mpc(1)*ngy + ngx
              ig_l = ig/npar_mpi
              ig_p = ig - ig_l*npar_mpi

              igrid_process(ngx, ngy, ngz) = ig_p
              igrid_l(ngx, ngy, ngz) = ig_l
              lgrid_l(ig_p) = lgrid_l(ig_p) + 1
           end do
        end do
     end do
  end if


  velo_solv_l(1:3, 1:inmpc%n_all_solv) = 0.0
  velo_mp_l(1:3, 1:nmp_real, irep) = 0.0


  n_solv = inmpc%n_all_solv
  ! -----------------------------------------------------------------
  !! av_xyz_coord_solv(1:3)=0.0
  !! av2_xyz_coord_solv(1:3)=0.0
  !! free-streaming process for mpc solvent particle
  tstep_colli = inmpc%nratio_colli_step*tstep
  do is = 1, n_solv
     xyz_solv_mpc(1:3, is) = xyz_solv_mpc(1:3, is) &
          + velo_solv_mpc(1:3, is)*tstep_colli
  end do
  

  ! -----------------------------------------------------------------
  ! grid (cell) devision process
  ! system particle and solvent particle of mpc is devided into cell
  TIME_S( tm_grid_mpc )
  call simu_grid_devision_mpc3(irep)
  TIME_E( tm_grid_mpc )

 
  ! -----------------------------------------------------------------
  ! rotate velocity
  TIME_S( tm_rotate_mpc )
!  call simu_rotate_velo_mpc2(velo_mp, irep)
  call simu_rotate_velo_mpc_rev3(velo_mp, irep, velo_solv_l, velo_mp_l)
  TIME_E( tm_rotate_mpc )


  ! -----------------------------------------------------------------
  ! rescale velocity
  TIME_S( tm_velo_mpc )
  call simu_velo_correct_simp_mpc3(velo_mp, irep, velo_solv_l, velo_mp_l)
  TIME_E( tm_velo_mpc )
  

end subroutine simu_collision_mpc3
