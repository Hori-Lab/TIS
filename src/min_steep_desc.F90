! min_steep_desc
!> @brief

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

subroutine min_steep_desc(flg_converge)

   use const_index
   use const_physical
   use const_maxsize

   use if_mloop
   use var_inp,    only : outfile
   use var_struct, only : nmp_all, xyz_mp_rep, pxyz_mp_rep
   use var_setp,   only : insimu
   use var_simu,   only : istep,istep_sim, force_mp, velo_mp, energy, energy_unit
   use var_emin,   only : inemin, lambda, norm_max,  lunout, func_check_lambda, func_check_norm
   use time, only : time_s, time_e, tm_energy, tm_force
   use mpiconst

   implicit none
   
   logical, intent(out) :: flg_converge  ! If converged, this flag become .TRUE.

   real(PREC) :: xyz(SDIM, nmp_all)
   real(PREC) :: pxyz(SDIM,nmp_all)
   real(PREC),save :: e_total
   real(PREC) :: e_total_new

   integer, parameter :: IREP = 1   ! Energy minimization is not available for REMD

   ! --------------------------------------------------------------
   ! Initialize (only at the 1st step)
   ! --------------------------------------------------------------
   if (istep == insimu%i_tstep_init) then
#ifdef MPI_PAR
      if (local_rank_mpi == 0) then
#endif
      if (inemin%i_out == 1) then
         lunout = outfile%data
      else if (inemin%i_out == 2) then
         lunout = outfile%opt
      endif
      if (inemin%i_out > 0) then
         write(lunout,'(a)',ADVANCE='no') '### Energy minimization : stage'
         write(lunout,'(i3)') istep_sim
         write(lunout,'(a)') '### The steepest descent method'
         write(lunout,'(a)') '###  istep ,   Energy     , Maximum norm , Lambda'
      endif
#ifdef MPI_PAR
      endif
#endif

      velo_mp(:,:,:) = 0.0e0_PREC

      TIME_S( tm_force )
      call force_sumup(force_mp, IREP)
      TIME_E( tm_force )

      ! Calc max(|Fn|)
      norm_max = max(maxval(force_mp), abs(minval(force_mp)))
      !Debug: write(*,*) 'norm_max=',norm_max

      TIME_S( tm_energy )
      call energy_sumup(IREP, velo_mp(:,:,IREP), energy(:,IREP), energy_unit(:,:,:,IREP))
      TIME_E( tm_energy )

      e_total = energy(E_TYPE%TOTAL, IREP)
      lambda = inemin%sd_lambda_init
      !Debug: write(*,*) 'step=1'
      !Debug: write(*,*) 'e_total=',e_total
   endif

   ! Save the current structure
   xyz(:,:) = xyz_mp_rep(:,:,IREP)
   pxyz(:,:) = pxyz_mp_rep(:,:,IREP)

   xyz_mp_rep(:,:,IREP) = xyz(:,:) + lambda * force_mp(:,:) / norm_max
   pxyz_mp_rep(:,:,IREP) = pxyz(:,:) + lambda * force_mp(:,:) / norm_max

   ! Calc energy
   TIME_S( tm_energy )
   call energy_sumup(IREP, velo_mp(:,:,IREP), energy(:,IREP), energy_unit(:,:,:,IREP))
   TIME_E( tm_energy )
   e_total_new = energy(E_TYPE%TOTAL, IREP)
   !Debug: write(*,*) 'e_total=',e_total
   !Debug: write(*,*) 'e_total_new=',e_total_new

   do while (e_total_new > e_total)
      lambda = lambda * inemin%sd_rho_reject

      ! Check the convergence
      if (func_check_lambda()) then
         flg_converge = .TRUE.
         xyz_mp_rep(:,:,IREP) = xyz(:,:)
         pxyz_mp_rep(:,:,IREP) = pxyz(:,:)
      endif

      ! Update to a candidate structure: 
      !         xyz(n+1) = xyz(n) + lambda * F(n) / norm_max
      !         where norm_max ||F(n)|| is max(|F(n)|)
      xyz_mp_rep(:,:,IREP) = xyz(:,:) + lambda * force_mp(:,:) / norm_max
      pxyz_mp_rep(:,:,IREP) = pxyz(:,:) + lambda * force_mp(:,:) / norm_max

      ! Calc energy
      TIME_S( tm_energy )
      call energy_sumup(IREP, velo_mp(:,:,IREP), energy(:,IREP), energy_unit(:,:,:,IREP))
      TIME_E( tm_energy )
      e_total_new = energy(E_TYPE%TOTAL, IREP)
      !Debug: write(*,*) 'Rejected: lambda=',lambda
   enddo

   ! Output 
#ifdef MPI_PAR
   if (local_rank_mpi == 0) then
#endif
   if (inemin%i_out > 0 .AND. mod(istep, inemin%n_out_step) == 0) then
      write(lunout, '(i10,1x,g15.8,1x,g15.8,1x,g15.8)') istep, e_total_new, norm_max, lambda
   endif
#ifdef MPI_PAR
   endif
#endif

   e_total = e_total_new
   lambda = lambda * inemin%sd_rho_accept
   !Debug: write(*,*) 'Accepted: lambda=',lambda

   ! Calc force
   TIME_S( tm_force )
   call force_sumup(force_mp, IREP)
   TIME_E( tm_force )

   ! Calc max(|Fn|)
   norm_max = max(maxval(force_mp), abs(minval(force_mp)))
   !Debug: write(*,*) 'norm_max=',norm_max

   ! Check the convergence
   if (func_check_norm()) then
      flg_converge = .TRUE.
      return
   endif

endsubroutine min_steep_desc
