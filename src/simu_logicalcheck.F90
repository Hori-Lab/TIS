! simu_logicalcheck
!> @brief Check consistency in input file

subroutine simu_logicalcheck()

  use const_maxsize
  use const_index
  use var_setp,    only : inmisc
  use var_io,     only : i_run_mode
  use var_setp,    only : inmisc
  use var_replica
#ifdef MPI_PAR
  use mpiconst
#endif
  integer :: error_code
  character(CARRAY_MSG_ERROR) :: error_message

  error_code     = 0

  !-------------------------------
  ! Checking for i_run_mode
  !-------------------------------
  select case (i_run_mode)
  ! Debug Mode, Check the consistence between force and energy
  case (RUN%CHECK_FORCE)

  ! Constant temperature simulation 
  case (RUN%CONST_TEMP)

  ! Simulated annealing
  case (RUN%SA)
!     if (i_type /= REPTYPE%VOID) error_code = 10 ! i_type should be VOID

  ! Auto-search of T_f
  case (RUN%SEARCH_TF)
!     if (i_type /= REPTYPE%VOID) error_code = 20  ! i_type should be VOID

  ! Energy calculation at single point
  case (RUN%ENERGY_CALC)

  ! Energy calculation for DCD trajectory
  case (RUN%ENERGY_DCD)

  ! Replica exchange
  case (RUN%REPLICA)
    
  case (RUN%FMAT)

  case (RUN%EMIN)

  case (RUN%WIDOM)

  case (RUN%GC)

  case default
     error_code = 90

  endselect


  !-------------------------------
  ! Checking for i_type
  !-------------------------------
!  select case (i_type)
!  ! Not replica
!  case (REPTYPE%VOID)
!     if (n_replica /= 1)             error_code = 110
!
!  ! Temperature replica
!  case (REPTYPE%TEMP)
!     if (i_run_mode /= RUN%REPLICA)   error_code = 210
!     if (n_replica <= 1)              error_code = 220
!     if (n_replica > MXREPLICA)       error_code = 222
!     if (.not. flg_rep%temp)          error_code = 230
!
!  ! Ion strength replica
!  case (REPTYPE%ION)
!     if (i_run_mode /= RUN%REPLICA)   error_code = 310
!     if (n_replica <= 1)              error_code = 320
!     if (n_replica > MXREPLICA)       error_code = 322
!     if (.not. flg_rep%ion)           error_code = 330
!
!  case (REPTYPE%TEMP_ION)
!     if (i_run_mode /= RUN%REPLICA)   error_code = 410
!     if (n_replica <= 1)              error_code = 420
!     if (n_replica > MXREPLICA)       error_code = 422
!     if (.not. flg_rep%temp)          error_code = 430
!     if (.not. flg_rep%ion)           error_code = 432
!
!  case default
!     error_code = 900
!
!  endselect

  ! i_reset_struct is only available either FMAT or REPLICA
  if (inmisc%i_reset_struct == 1) then
     if ((i_run_mode == RUN%FMAT)      .OR.&
         (i_run_mode == RUN%SEARCH_TF) .OR.&
        (i_run_mode == RUN%REPLICA .AND. inrep%flg_opt_temp)) then
        continue
     else
        write(error_message, *) 'i_reset_struct is not available with current i_run_mode'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  endif

  !-------------------------------
  ! output error message, then program stop
  if (error_code /= 0) then
     write(error_message,*) 'error_code = ',error_code
     call util_error(ERROR%WARN_ALL,error_message)
     write(error_message,*) 'logical error, PROGRAM STOP'
     call util_error(ERROR%STOP_ALL, error_message)
  endif
  
endsubroutine simu_logicalcheck
