! min_conj_grad
!> @brief  This subroutine performes a energy minimization by the conjugated gradient method
!          (to be more specific, Polak-Ribiere(PR) method), with a line search approach.
!
! Reference:
!    Nocedal J and Wright SJ, "Numerical Optimization", Springer-Verlag (New York), 1999

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

subroutine min_conj_grad(flg_converge)

   use const_index
   use const_physical
   use const_maxsize

   use if_mloop
   use var_io,    only : outfile
   use var_struct, only : nmp_all, xyz_mp_rep, pxyz_mp_rep
   use var_setp,   only : insimu
   use var_simu,   only : istep,istep_sim, force_mp, velo_mp, energy, energy_unit, dxyz_mp
   use var_emin,   only : inemin, lambda, norm_max, lunout, pvec, func_check_lambda, func_check_norm
   use time, only : time_s, time_e, tm_neighbor, tm_energy, tm_force
   use mpiconst

   implicit none
   
   logical, intent(out) :: flg_converge  ! If converged, this flag become .TRUE.

   real(PREC) :: beta
   real(PREC) :: dotp1, dotp2
   real(PREC) :: xyz(SDIM, nmp_all)
   real(PREC) :: pxyz(SDIM,nmp_all)
   real(PREC) :: dxyz(SDIM,nmp_all)
   real(PREC) :: force_new(SDIM,nmp_all)
   real(PREC) :: e_total_new
   real(PREC), save :: e_total
   integer, parameter :: IREP = 1   ! Energy minimization is not available for REMD

   lambda = inemin%cg_lambda_init

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
         write(lunout,'(a)') '### The conjugate gradient method'
         write(lunout,'(a)') '###  istep ,   Energy     , Maximum norm , Lambda'
      endif
#ifdef MPI_PAR
      endif
#endif

      velo_mp(:,:,:) = 0.0e0_PREC

      TIME_S( tm_energy )
      call energy_sumup(IREP, velo_mp(:,:,IREP), energy(:,IREP), energy_unit(:,:,:,IREP))
      TIME_E( tm_energy )
      e_total = energy(E_TYPE%TOTAL, IREP)

      TIME_S( tm_force )
      call force_sumup(force_mp(:,:,IREP), IREP)
      TIME_E( tm_force )

      allocate(pvec(SDIM, nmp_all))
      pvec(:,:) = force_mp(:,:,IREP)
      norm_max = max(maxval(pvec), abs(minval(pvec)))
   endif
 
   ! Save the current structure
   xyz(:,:) = xyz_mp_rep(:,:,IREP)
   pxyz(:,:) = pxyz_mp_rep(:,:,IREP)
   dxyz(:,:) = dxyz_mp(:,:,IREP)

   ! Update to a candidate structure: 
   xyz_mp_rep(:,:,IREP) = xyz(:,:) + lambda * pvec(:,:)
   pxyz_mp_rep(:,:,IREP) = pxyz(:,:) + lambda * pvec(:,:)
   dxyz_mp(:,:,IREP) = dxyz(:,:) + lambda * pvec(:,:)

   ! Calc energy
   TIME_S( tm_energy )
   call energy_sumup(IREP, velo_mp(:,:,IREP), energy(:,IREP), energy_unit(:,:,:,IREP))
   TIME_E( tm_energy )
   e_total_new = energy(E_TYPE%TOTAL, IREP)

   ! Calc force
   TIME_S( tm_force )
   call force_sumup(force_new, IREP)
   TIME_E( tm_force )

   dotp1 = func_dotp(-force_mp(:,:,IREP), pvec)
   dotp2 = func_dotp(-force_new, pvec)

   ! "Strong Wolfe" condition is tested.
   do while ((e_total_new > e_total + inemin%cg_wolfe_c1 * lambda * dotp1) .OR. &
             (abs(dotp2) > inemin%cg_wolfe_c2 * abs(dotp1)               )  )
      ! Decrease lambda
      lambda = inemin%cg_rho * lambda

      ! Check the convergence
      if (func_check_lambda()) then
         flg_converge = .TRUE.
         xyz_mp_rep(:,:,IREP) = xyz(:,:)
         pxyz_mp_rep(:,:,IREP) = pxyz(:,:)
         dxyz_mp(:,:,IREP) = dxyz(:,:)
         return
      endif
      ! Update to a candidate structure: 
      xyz_mp_rep(:,:,IREP) = xyz(:,:) + lambda * pvec(:,:)
      pxyz_mp_rep(:,:,IREP) = pxyz(:,:) + lambda * pvec(:,:)
      dxyz_mp(:,:,IREP) = dxyz(:,:) + lambda * pvec(:,:)
   
      ! Calc energy
      TIME_S( tm_energy )
      call energy_sumup(IREP, velo_mp(:,:,IREP), energy(:,IREP), energy_unit(:,:,:,IREP))
      TIME_E( tm_energy )
      e_total_new = energy(E_TYPE%TOTAL, IREP)
   
      ! Calc force
      TIME_S( tm_force )
      call force_sumup(force_new, IREP)
      TIME_E( tm_force )
   
      !dotp1 = func_dotp(-force_mp, pvec)
      dotp2 = func_dotp(-force_new, pvec)
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

   ! FR (Fletcher-Reeves) method
   !dotp1 = func_dotp(force_new, force_new)
   !dotp2 = func_dotp(force_mp, force_mp)

   ! PR (Polak-Ribiere) method
   dotp1 = func_dotp(-force_new, (force_mp(:,:,IREP) - force_new))
   dotp2 = func_dotp(force_mp(:,:,IREP), force_mp(:,:,IREP))

   ! modification for PR method
   beta = max(dotp1 / dotp2, 0.0e0_PREC)

   pvec(:,:) = force_new(:,:) + beta * pvec(:,:)
   norm_max = max(maxval(pvec), abs(minval(pvec)))

   ! Check the convergence
   if (func_check_norm()) then
      flg_converge = .TRUE.
      return
   endif

   e_total = e_total_new
   force_mp(:,:,IREP) = force_new(:,:)


! ###########################################################################
contains
   
   ! Calculate dot product of "a" and "b",
   ! assuming "a" and "b" are one-dimensional vector although they are two-dimensional array.
   real(PREC) function func_dotp(a,b)
      implicit none
      real(PREC), intent(in) :: a(:,:)
      real(PREC), intent(in) :: b(:,:)
      integer :: i,j
      ! CAUTION: 
      ! This function assumes that the size (and shape) of "a" and "b" are identical.

      func_dotp = 0.0e0_PREC
      do j = lbound(a,2), ubound(a,2)
         do i = lbound(a,1), ubound(a,1)
            func_dotp = func_dotp + a(i,j) * b(i,j)
         enddo
      enddo
   endfunction func_dotp

endsubroutine min_conj_grad
