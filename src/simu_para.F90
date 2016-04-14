subroutine simu_para()

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : i_run_mode
  use var_setp,    only : insimu, inmmc
  use var_replica, only : inrep, n_replica_all, lab2step, step_ratio
  use var_struct,  only : nmp_real, cmass_mp
  use var_simu,    only : istep_sim, istep, mstep, ntstep, tstep, &
                          tempk, nstep_opt_temp, ntstep_max, &
                          rcmass_mp, tstep2, tsteph, &
                          qscore, n_exchange, max_exchange, &
                          em_mid, em_depth, em_sigma
  use var_fmat,    only : infmat, fmat_clear

  implicit none

  ! -----------------------------------------------------------------
  ! local variables
  integer :: irep
  character(CARRAY_MSG_ERROR)   :: error_message


  ! -----------------------------------------------------------------
  ! setting parameters
  tstep  = insimu%tstep_size
  tempk  = insimu%tempk
  nstep_opt_temp = inrep%n_step_opt_temp

  rcmass_mp(1:nmp_real) = 1.0e0_PREC / cmass_mp(1:nmp_real)
  tstep2           = 0.5e0_PREC * tstep * tstep
  tsteph           = 0.5e0_PREC * tstep

  ! mcanonical  : setting parameters
  if(inmmc%i_modified_muca == 1)then
     em_depth = inmmc%em_depth
     em_mid = inmmc%em_mid
     em_sigma = inmmc%em_sigma
  endif


  mstep = 1  ! for not searching Tf

  ! anealing or heating
  if(i_run_mode == RUN%SA) then
     ntstep = insimu%n_tstep(istep_sim)

     ! # of replica should be one
     if (n_replica_all /= 1) then
        write(error_message,*) 'defect at simu_para, PROGRAM STOP'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

     istep = 0
     call simu_anneal(istep, tempk)

  ! searchingtf
  else if(i_run_mode == RUN%SEARCH_TF) then
     ntstep = insimu%n_tstep(istep_sim)

     ! # of replica should be one
     if (n_replica_all /= 1) then
        write(error_message,*) 'defect at simu_para, PROGRAM STOP'
        call util_error(ERROR%STOP_ALL, error_message) 
     endif

     istep = 0
     call simu_searchingtf(istep, ntstep, qscore, tempk)
     mstep = MXSEARCHINGTF    ! defined by const_maxsize.F90 (default=1000)

  ! fluctuation matching
  else if (i_run_mode == RUN%FMAT) then
     call fmat_clear()
     ntstep = infmat%n_step

  ! Replica exchange method
  else if (i_run_mode == RUN%REPLICA) then
     ntstep = insimu%n_tstep(istep_sim)

     if (inrep%flg_opt_temp) then
        ! When the number of stages is specifed for FO-REM
        if (inrep%n_stage_opt_temp /= MX_REP_OPT_STAGE) then
           ntstep = MX_NTSTEP
        endif
     endif

     if (inrep%i_loadbalance >= 1) then
        max_exchange = 0
        do irep = 1, n_replica_all
          max_exchange = max(max_exchange, ntstep / lab2step(irep))
        enddo
        ntstep_max = max_exchange * maxval(lab2step) * step_ratio
        ntstep = ntstep_max
        write(6,*) ' - maxval(lab2step) - ', maxval(lab2step)
        write(6,*) ' - ntstep, max_exchange, ntstep_max - ', ntstep, max_exchange, ntstep_max
     endif

     n_exchange = 0

  else
     ntstep = insimu%n_tstep(istep_sim)

  endif

  ! check the limitation of total step number
  if (ntstep > MX_NTSTEP) then
     ntstep = MX_NTSTEP
     write(error_message,'(a,x,i25)') &
        'the # of total step is over the limit. "ntstep" is decreased to',ntstep
     call util_error(ERROR%WARN_ALL, error_message)
  endif

end subroutine simu_para
