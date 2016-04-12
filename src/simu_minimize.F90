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
   use var_setp, only : insimu
   use var_simu, only : istep,istep_sim
   use var_emin, only : inemin

   use time, only : time_s, time_e, tm_neighbor

   implicit none
   
   logical, intent(out) :: flg_converge  ! When converged, this flag become .TRUE.
   integer, parameter :: IREP = 1   ! dummy: Energy minimization is not available for REMD
   character(CARRAY_MSG_ERROR) :: error_message

   ! --------------------------------------------------------------
   ! Calc neighbour list
   ! --------------------------------------------------------------
   if(mod(istep, insimu%n_step_neighbor) == 1 .OR. istep == insimu%i_tstep_init) then  
      TIME_S( tm_neighbor )
      call simu_neighbor(IREP)
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
