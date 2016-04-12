! ***********************************************************************
subroutine simu_velo_correct_simp_mpc3(velo_mp, irep, velo_solv_l, velo_mp_l)

  use const_maxsize, only : PREC, MXMP, MXSOLV_MPC
  use const_physical, only : BOLTZC
  use var_setp, only : ifix_mp
  use var_struct, only : nmp_real, cmass_mp, xyz_mp_rep
  use var_simu, only : istep, tempk
  use var_mpc, only : inmpc, cmass_solv_mpc, velo_solv_mpc, xyz_solv_mpc, &
       igrid_process, lgrid_l, igrid_l, &
       nis_l, nimp_l, isolv2grid_mpc_l, imp2grid_mpc_l

  use var_replica, only : n_replica_all, n_replica_mpi

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none
  
  ! --------------------------------------------------------------------
  real(PREC), intent(inout) :: velo_mp(:,:,:)
  real(PREC), intent(inout) :: velo_solv_l(:, :), velo_mp_l(:, :, :)
  integer,    intent(in) :: irep

  ! --------------------------------------------------------------------
  !local variables
  integer :: is, imp, n_solv
  !integer :: is_l, imp_l
  real(PREC) :: P_total_mpc(3), velo_cm_mpc(3), total_mass_mpc
  real(PREC) :: rescale_fact_All, total_kine_E2_All
  integer :: i_count_mp_All

  ! --------------------------------------------------------------------
#ifdef MPI_PAR
  call mpi_allreduce(velo_solv_l, velo_solv_mpc, 3*MXSOLV_MPC, &
       PREC_MPI, MPI_SUM, mpi_comm_local, ierr)
  call mpi_allreduce(velo_mp_l, velo_mp, 3*nmp_real*n_replica_mpi, &
       PREC_MPI, MPI_SUM, mpi_comm_local, ierr)
!  do imp = 1, 10
!     write (*, *) myrank, imp, velo_mp_l(1, imp, irep), velo_mp(1, imp, irep)
!  end do
!  stop
#else
  velo_solv_mpc(1:3, 1:inmpc%n_all_solv) = velo_solv_l(1:3, 1:inmpc%n_all_solv)
  velo_mp(1:3, 1:nmp_real, irep) = velo_mp_l(1:3, 1:nmp_real, irep)
#endif

  ! --------------------------------------------------------------------
  !!------------------------------------------
  !! total-linear momentum should be zero.
  !!------------------------------------------
  total_mass_mpc = 0.0e0_PREC
  P_total_mpc(1:3) = 0.0e0_PREC

  n_solv = inmpc%n_all_solv
  do is = 1, n_solv
!  do is_l = 1, nis_l
!     is = isolv2grid_mpc_l(2, is_l)
     P_total_mpc(1:3) = P_total_mpc(1:3) &
          + cmass_solv_mpc(is) * velo_solv_mpc(1:3, is)
     total_mass_mpc = total_mass_mpc + cmass_solv_mpc(is)
  end do
  
  do imp = 1, nmp_real
!  do imp_l = 1, nimp_l
!     imp = imp2grid_mpc_l(2, imp_l)
     if(ifix_mp(imp) == 1) cycle
     P_total_mpc(1:3) = P_total_mpc(1:3) &
          + cmass_mp(imp) * velo_mp(1:3, imp, irep)
     total_mass_mpc = total_mass_mpc + cmass_mp(imp)
  end do
  
  velo_cm_mpc(1:3) = P_total_mpc(1:3) / total_mass_mpc
  
  !!velocity correction
  do is = 1, n_solv
!  do is_l = 1, nis_l
!     is = isolv2grid_mpc_l(2, is_l)
     velo_solv_mpc(1:3, is) = velo_solv_mpc(1:3, is) - velo_cm_mpc(1:3)
  end do

  do imp = 1, nmp_real
!  do imp_l = 1, nimp_l
!     imp = imp2grid_mpc_l(2, imp_l)
     if(ifix_mp(imp) == 1) cycle
     velo_mp(1:3, imp, irep) = velo_mp(1:3, imp, irep) - velo_cm_mpc(1:3)
  end do
  
     
  if(mod(istep, inmpc%nratio_vcorrect_step) == 0) then

     !!------------------------------------------
     !!calculate rescaling factor for solvent particle velocity to reset T
     !!-----------------------------------------
     total_kine_E2_All = 0.0e0_PREC
     i_count_mp_All = 0
        
     !!<<for mpc solvent particle>>
     do is = 1, n_solv
!     do is_l = 1, nis_l
!        is = isolv2grid_mpc_l(2, is_l)
        total_kine_E2_All = total_kine_E2_All + cmass_solv_mpc(is) &
             * (velo_solv_mpc(1, is)**2 &
             + velo_solv_mpc(2, is)**2 + velo_solv_mpc(3, is)**2)
     end do
     i_count_mp_All = i_count_mp_All + 3*n_solv

     rescale_fact_All = sqrt((BOLTZC*tempk)*(i_count_mp_All - 3)/total_kine_E2_All)
     
     !!-----------------------------------------
     !! velocity correct
     !!-----------------------------------------
     do is = 1, n_solv
!     do is_l = 1, nis_l
!        is = isolv2grid_mpc_l(2, is_l)
        velo_solv_mpc(1:3, is) = velo_solv_mpc(1:3, is) * rescale_fact_All
     end do
  end if

end subroutine simu_velo_correct_simp_mpc3
