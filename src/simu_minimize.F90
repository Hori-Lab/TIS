! simu_minimize
!> @brief

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

subroutine simu_minimize(flg_converge)

   use const_index
   use const_maxsize

   use if_mloop
   use var_setp, only : insimu, inmisc, inpara
   use var_simu, only : istep,istep_sim, dxyz_mp
   use var_struct,only: nmp_real
   use var_emin, only : inemin
   use var_io,   only : ifile_out_neigh, outfile
   use time, only : time_s, time_e, tm_neighbor

   implicit none
   
   logical, intent(out) :: flg_converge  ! When converged, this flag become .TRUE.
   integer, parameter :: IREP = 1   ! dummy: Energy minimization is not available for REMD
   character(CARRAY_MSG_ERROR) :: error_message
   integer :: imp
   real(PREC) :: d2, d2max, d2max_2nd

   ! --------------------------------------------------------------
   ! Calc neighbour list
   ! --------------------------------------------------------------
   if (inmisc%i_neigh_dynamic == 1) then
      TIME_S( tm_neighbor )
 
      d2max = 0.0e0_PREC
      d2max_2nd = 0.0e0_PREC
      do imp = 1, nmp_real
         d2 = dot_product(dxyz_mp(1:3,imp,irep), dxyz_mp(1:3,imp,irep))
         if (d2 > d2max) then
            d2max_2nd = d2max
            d2max = d2
         else if (d2 > d2max_2nd) then
            d2max_2nd = d2
         else
            cycle
         endif
      enddo

      if (sqrt(d2max) + sqrt(d2max_2nd) > inpara%neigh_margin) then
         if (ifile_out_neigh == 1) then
            write(outfile%neigh, '(i10,1x,i5,1x,f4.1,1x,f4.1,1x,f4.1)',advance='no') &
                                 istep, irep, d2max, d2max_2nd, d2max+d2max_2nd
         endif
         call neighbor(irep)
         dxyz_mp(:,:,irep) = 0.0e0_PREC
      endif
 
      TIME_E( tm_neighbor )
   else if(mod(istep, insimu%n_step_neighbor) == 1 .OR. istep == insimu%i_tstep_init) then  
      TIME_S( tm_neighbor )
      if (ifile_out_neigh == 1) then
         write(outfile%neigh, '(i10,1x,i5)',advance='no') istep, irep
      endif
      call neighbor(IREP)
      TIME_E( tm_neighbor )
   end if
 
   ! --------------------------------------------------------------
   ! Main procedure
   ! --------------------------------------------------------------
   select case (inemin%i_method(istep_sim))

   ! Steepest Descent
   case(EMIN_METHOD%SD)
      call min_steep_desc(flg_converge)

   ! Conjugate Gradient
   case(EMIN_METHOD%CG)
      call min_conj_grad(flg_converge)

   case default
      error_message = 'Error: unknown method for energy minimization'
      call util_error(ERROR%STOP_ALL, error_message)
 
   endselect
 
endsubroutine simu_minimize
