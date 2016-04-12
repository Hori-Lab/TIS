!simu_collision_mpc2
!> @brief

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

subroutine simu_collision_mpc2(irep)

  use const_maxsize
  use const_physical
!  use const_index
  use if_mloop
!  use if_write
!  use if_energy
  use var_struct, only : nmp_real, xyz_mp_rep, cmass_mp
  use var_simu, only : istep, velo_mp, tempk
  use var_replica, only : n_replica_mpi
  use var_mpc, only : inmpc

  use time, only : tm_grid_mpc, tm_rotate_mpc, tm_velo_mpc, &
                   time_s, time_e

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! -----------------------------------------------------------------
  integer, intent(in) :: irep


  ! -----------------------------------------------------------------
  ! local variables
!  integer :: is
!  real(PREC) :: tstep_colli

  ! -----------------------------------------------------------------
  if(mod(istep, inmpc%nratio_colli_step) /= 0) then
     return
  end if
  
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  TIME_S( tm_velo_mpc )
  call simu_stream_mpc2(velo_mp, irep)
  TIME_E( tm_velo_mpc )

  ! -----------------------------------------------------------------
  !! av_xyz_coord_solv(1:3)=0.0
  !! av2_xyz_coord_solv(1:3)=0.0
  !! free-streaming process for mpc solvent particle
!  n_solv = inmpc%n_all_solv
!  tstep_colli = inmpc%nratio_colli_step*tstep
!  do is = 1, n_solv
!     xyz_solv_mpc(1:3, is) = xyz_solv_mpc(1:3, is) &
!          + velo_solv_mpc(1:3, is)*tstep_colli
!  end do
  

  ! -----------------------------------------------------------------
  ! grid (cell) devision process
  ! system particle and solvent particle of mpc is devided into cell
  TIME_S( tm_grid_mpc )
  call simu_grid_devision_mpc2(irep)
  TIME_E( tm_grid_mpc )
 

  ! -----------------------------------------------------------------
  ! rotate velocity
  TIME_S( tm_rotate_mpc )
!  call simu_rotate_velo_mpc2(velo_mp, irep)
  call simu_rotate_velo_mpc_rev2(velo_mp, irep)
  TIME_E( tm_rotate_mpc )


  ! -----------------------------------------------------------------
  ! rescale velocity
!  TIME_S( tm_velo_mpc )
!  call simu_velo_correct_simp_mpc2(velo_mp, irep)
!  TIME_E( tm_velo_mpc )
  
#ifdef MPI_PAR
  end if

  call MPI_Bcast (velo_mp, SPACE_DIM*nmp_real*n_replica_mpi, PREC_MPI,0,MPI_COMM_WORLD,ierr)

!  call MPI_Bcast (velo_solv_mpc, 3*MXSOLV_MPC, PREC_MPI,0,MPI_COMM_WORLD,ierr)
#endif


end subroutine simu_collision_mpc2
