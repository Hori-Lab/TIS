! inp_replica_para
!> @brief Read parameters from .inp file for doing replica exchange method

subroutine inp_replica_para()

  use const_maxsize
  use const_index
  use const_physical
  use var_io,     only : infile, i_run_mode
  use var_replica, only : inrep, flg_rep, flg_npar_rep, &
                          n_replica_all, n_dimension
  use var_setp, only : inmisc

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  integer :: icol
  integer :: irep, ivar, i_exchange, i_opt_temp, i_pull
  integer :: imp, jmp
  real(PREC) :: coef, length

  integer :: igrp
  real(PREC) :: z, kxy, kz

  character(4)  :: kfind      
  character(16) :: char_query
  character(CARRAY_MXCOLM)  :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM)  :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MXCOLM) :: ctmp00, ctmp02
  character(CARRAY_MXCOLM) :: cdummy
  character(CARRAY_MSG_ERROR) :: error_message

#ifdef _DEBUG
  write(6,*) 'inp_replica_para : START'
#endif
   
  luninp = infile%inp
  lunout  = 6   ! At this moment, '.data' file is not exist.
  
  ! default
  i_exchange = -1
  i_opt_temp = 0
  n_replica_all  = 1
  inrep%npar_rep = 1
  n_dimension    = 0

  ! initialize
  flg_npar_rep = .false.
  inrep%flg_opt_temp = .false.
  inrep%n_step_opt_temp = -1
  inrep%n_stage_opt_temp = MX_REP_OPT_STAGE
  inrep%n_step_replica = -1
  inrep%n_step_save    = -1
  inrep%n_replica(1:REPTYPE%MAX) = -1
  inrep%i_style(1:REPTYPE%MAX)  = REPVARSTYLE%VOID
!  inrep%i_loadbalance  = -1
!  inrep%n_step_adjust  = -1
!  inrep%n_adjust_interval  = -1
  inrep%lowest(1:REPTYPE%MAX)   = INVALID_VALUE
  inrep%highest(1:REPTYPE%MAX)  = INVALID_VALUE
  !inrep%interval(1:REPTYPE%MAX) = INVALID_VALUE
  inrep%var(1:MXREPLICA, 1:REPTYPE%MAX) = INVALID_VALUE
  inrep%flg_exchange = .true.
  inrep%exponent_alpha(1:REPTYPE%MAX) = 1.0
  inrep%n_period_prob = 2
  flg_rep(1:REPTYPE%MAX) = .false.

  if (i_run_mode /= RUN%REPLICA) then
     inrep%n_replica(1:REPTYPE%MAX) = 1
     inrep%n_step_replica = 0
     inrep%n_step_save    = 0  ! default
     return
  end if

  ! Pulling
  i_pull = 0
  inrep%pull_direction(1:SDIM, 1:MXPULLING) = 0.0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  !######################################################################
  rewind(luninp)
     
  !-----------------------------------------------
  ! Reading "replica" field
  !-----------------------------------------------
  call ukoto_uiread2(luninp, lunout, 'replica         ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "replica" field in the input file, PROGRAM STOP'
     call util_error(ERROR%STOP_ALL, error_message)
  endif
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     do iequa = 1, nequat
        cvalue = 'n_replica_temp'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inrep%n_replica(REPTYPE%TEMP), cvalue)
        
        cvalue = 'i_exchange'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             i_exchange, cvalue)
        
        cvalue = 'i_opt_temp'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             i_opt_temp, cvalue)
        
        cvalue = 'n_step_opt_temp'
        call ukoto_lvalue2(lunout, csides(1, iequa), &
             inrep%n_step_opt_temp, cvalue)
        
        cvalue = 'n_stage_opt_temp'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inrep%n_stage_opt_temp, cvalue)
        
        cvalue = 'n_replica_ion'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inrep%n_replica(REPTYPE%ION), cvalue)
        
        cvalue = 'n_replica_pull'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inrep%n_replica(REPTYPE%PULL), cvalue)

        cvalue = 'n_replica_wind'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inrep%n_replica(REPTYPE%WIND), cvalue)

!        cvalue = 'n_replica_winz'
!        call ukoto_ivalue2(lunout, csides(1, iequa), &
!             inrep%n_replica(REPTYPE%WINZ), cvalue)
        
        cvalue = 'n_step_exchange'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inrep%n_step_replica, cvalue)
        
        cvalue = 'n_step_save_rep'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inrep%n_step_save, cvalue)
        
        cvalue = 'npar_rep'
        if(csides(1, iequa) == 'npar_rep') then
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inrep%npar_rep, cvalue)
           flg_npar_rep = .true.
        end if
     end do
  end do
  
  
#ifdef MPI_PAR
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     do iequa = 1, nequat
        
!        cvalue = 'i_loadbalance'
!        call ukoto_ivalue2(lunout, csides(1, iequa), &
!             inrep%i_loadbalance, cvalue)
        
!        cvalue = 'n_step_adjust'
!        call ukoto_ivalue2(lunout, csides(1, iequa), &
!             inrep%n_step_adjust, cvalue)
!        
!        cvalue = 'n_adjust_interval'
!        call ukoto_ivalue2(lunout, csides(1, iequa), &
!             inrep%n_adjust_interval, cvalue)
        
        cvalue = 'n_period_prob'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inrep%n_period_prob, cvalue)
     end do
  end do
#endif
  
!  ! checking input variables
!  if(inrep%i_loadbalance /= -1) then
!     if(inrep%npar_rep == 1) then
!        inrep%i_loadbalance = 0
!     end if
!  else
!     inrep%i_loadbalance = 0
!  end if
  
  if (i_exchange == -1) then
     inrep%flg_exchange = .true.  ! default (i_exchange is not specified in the input)
  else if (i_exchange == 1) then
     inrep%flg_exchange = .true.
  else if (i_exchange == 0) then
     inrep%flg_exchange = .false.
  else
     error_message = 'Error: in reading i_exchage, PROGRAM STOP'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  if (i_opt_temp == 1) then
     inrep%flg_opt_temp = .true.
  else if (i_opt_temp == 0) then
     inrep%flg_opt_temp = .false.
  else
     error_message = 'Error: in reading i_opt_temp, PROGRAM STOP'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  ! n_step_replica
  if(inrep%n_step_replica < 1) then
     error_message = 'Error: invalid value for replica:n_step_replica, PROGRAM STOP'
     call util_error(ERROR%STOP_ALL, error_message)
  endif
  ! n_step_save_rep
  if(inrep%n_step_save < 1) then
     inrep%n_step_save = inrep%n_step_replica  ! this is default
  endif

  ! set i_window flag which decide whether window force and energy will be calculated
  if (inrep%n_replica(REPTYPE%WIND) > 0) then
     inmisc%i_window = 1
  else
     inmisc%i_window = 0
  endif

!  ! set i_window flag which decide whether window force and energy will be calculated
!  if (inrep%n_replica(REPTYPE%WINZ) > 0) then
!     inmisc%i_winz = 1
!  else
!     inmisc%i_winz = 0
!  endif
  
  ! set flg_rep, n_dimension, and n_replica_all
  do ivar = 1, REPTYPE%MAX
     if (inrep%n_replica(ivar) <= 0) then
        inrep%n_replica(ivar) = 1
        cycle
        
     elseif (inrep%n_replica(ivar) > MXREPLICA) then
        write(error_message,*)'Error: invalid value for replica: n_replica_',CHAR_REPTYPE(ivar),', PROGRAM STOP'
        call util_error(ERROR%STOP_ALL, error_message)
        
     else
        flg_rep(ivar) = .true.
        
        n_dimension = n_dimension + 1
        if (n_dimension > MXREPDIM) then
           error_message = 'Error: invalid value for replica: too many dimensions, PROGRAM STOP'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        
        n_replica_all = n_replica_all * inrep%n_replica(ivar) 
        if (n_replica_all > MXREPLICA) then
           error_message = 'Error: invalid value for replica: (#replica) > MXREPLICA, PROGRAM STOP'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        
     endif
  enddo

  if (flg_rep(REPTYPE%PULL)) then
     inmisc%i_pulling = 1
  else
     inmisc%i_pulling = 0
  endif
  
  write(lunout, *) 'i_exchange = ', inrep%flg_exchange
  write(lunout, *) 'n_replica_temp  = ', inrep%n_replica(REPTYPE%TEMP)
  write(lunout, *) 'n_replica_ion   = ', inrep%n_replica(REPTYPE%ION)
  write(lunout, *) 'n_replica_pull  = ', inrep%n_replica(REPTYPE%PULL)
  write(lunout, *) 'n_replica_wind  = ', inrep%n_replica(REPTYPE%WIND)
  write(lunout, *) 'n_step_exchange = ', inrep%n_step_replica
  write(lunout, *) 'n_step_save_rep = ', inrep%n_step_save
  write(lunout, *) 'i_opt_temp      = ', inrep%flg_opt_temp
  write(lunout, *) 'n_step_opt_temp = ', inrep%n_step_opt_temp
  write(lunout, *) 'n_stage_opt_temp= ', inrep%n_stage_opt_temp
  write(lunout, *) 'npar_rep = ', inrep%npar_rep
  
  if (n_dimension < 1) then
     error_message = 'Error: invalid value for replica: #dimension < 1, PROGRAM STOP'
     call util_error(ERROR%STOP_ALL, error_message)
  endif
  
  if (inrep%flg_opt_temp) then
     if (inrep%n_step_opt_temp <=0 ) then 
        error_message = 'Error: invalid value for n_step_opt_temp, PROGRAM STOP'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  endif
  
!  if (inrep%i_loadbalance >= 1) then
!     if (inrep%n_adjust_interval == -1) then
!        inrep%n_adjust_interval = 10    !  this is default
!        error_message = 'Error: invalid value for n_adjust_interval, PROGRAM STOP'
!        call util_error(ERROR%STOP_ALL, error_message)
!     endif
!  endif
  
  ! The explanation of replica parameter
!  write (lunout, *) 'i_loadbalance = ', inrep%i_loadbalance
!  if (inrep%i_loadbalance == 0) then
!     write (lunout, *) ' Without loadbalance(default) '
!  else if (inrep%i_loadbalance == 1) then
!     write (lunout, *) ' balanced by neighborlist_ele '
!  else if (inrep%i_loadbalance == 2) then
!     write (lunout, *) ' balanced by dynamic'
!  else
!     error_message = 'Error: invalid value about i_loadbalance'
!     call util_error(ERROR%STOP_ALL, error_message)
!  endif
!  if (inrep%i_loadbalance >=1) then
!     write (lunout, *) 'time step adjustment interval(x Replica exchange interval)  = ', &
!          inrep%n_adjust_interval
!  end if
  
  !-----------------------------------------------
  ! Reading "replica_XXX" field
  !-----------------------------------------------
  inrep%var(1:MXREPLICA,1:REPTYPE%MAX) = INVALID_VALUE
  
  do ivar = 1, REPTYPE%MAX
     
     if (.not. flg_rep(ivar)) then
        cycle
     endif
     
     rewind(luninp)
     
     write(char_query,'(a8,a4,a4)') 'replica_',CHAR_REPTYPE(ivar),'    '
     
     call ukoto_uiread2(luninp, lunout, char_query, kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if (kfind /= 'FIND') then
        error_message = 'Error: cannot find "'//char_query//'" field in the input file'
        !          error_message = 'Error: cannot find "'//char_query
        !          error_message = error_message//'" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
     
     do iline = 1, nlines
        ctmp00 = cwkinp(iline)
        
        call ukoto_uiequa2(lunout, ctmp00, nequat, csides)
        
        do iequa = 1, nequat
           ctmp02 = csides(1, iequa)
                      
           cvalue = 'i_style'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inrep%i_style(ivar), cvalue)
           
           if(ctmp02(1:7) == 'REPLICA') then
              irep = 0
              do icol = 9, CARRAY_MXCOLM
                 if(ctmp02(icol:icol) == ')') exit
              end do
              read (ctmp02(9:icol-1), *) irep
              if (irep < 1 .OR. irep > inrep%n_replica(ivar)) then
                 write(error_message,*) 'Error: in reading REPLICA(', irep, ')', inrep%n_replica(ivar), 'PROGRAM STOP'
                 call util_error(ERROR%STOP_ALL, error_message)
              endif
              
              cvalue = ctmp02(1:icol)
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   inrep%var(irep, ivar), cvalue)
           end if

           cvalue = 'value_lowest'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inrep%lowest(ivar), cvalue)
           
           cvalue = 'value_highest'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inrep%highest(ivar), cvalue)
           
           cvalue = 'exponent'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inrep%exponent_alpha(ivar), cvalue)
        end do ! iequa

        if (ctmp00(1:4) == 'WINZ') then
           irep = 0

           do icol = 6, CARRAY_MXCOLM
              if (ctmp00(icol:icol) == ')') exit
           end do

           read (ctmp00(6:icol-1), *) irep

           if (irep < 1 .OR. irep > inrep%n_replica(ivar)) then
              write(error_message,*) 'Error: in reading WINZ(', irep, ')', inrep%n_replica(ivar), 'PROGRAM STOP'
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           write(*, *) ctmp00(icol+1:CARRAY_MXCOLM)
           
           read(ctmp00(icol+1:CARRAY_MXCOLM), *) igrp, z, kxy, kz

!           inrep%winz_igrp(irep) = igrp
!           inrep%winz_z(irep) = z
!           inrep%winz_kxy(irep) = kxy
!           inrep%winz_kz(irep) = kz

           inrep%var(irep, ivar) = irep

        end if
           
        if (ctmp00(1:6) == 'WINDOW') then
           irep = 0
           
           do icol = 8, CARRAY_MXCOLM
              if (ctmp00(icol:icol) == ')') exit
           end do
           
           read (ctmp00(8:icol-1), *) irep
           
           if (irep < 1 .OR. irep > inrep%n_replica(ivar)) then
              write(error_message,*) 'Error: in reading WINDOW(', irep, ')', inrep%n_replica(ivar), 'PROGRAM STOP'
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           
           read(ctmp00(icol+1:CARRAY_MXCOLM), *) imp, jmp ,coef, length
           
           inrep%window_property(irep, WINDTYPE%COEF) = coef
           inrep%window_property(irep, WINDTYPE%LENGTH) = length
           inrep%window_mp_id(irep, WINDTYPE%IMP) = imp
           inrep%window_mp_id(irep, WINDTYPE%JMP) = jmp
           
           inrep%var(irep, ivar) = irep
        end if

        if (ivar==REPTYPE%PULL .and. ctmp00(1:15)=='PULL_UNRAVEL_CF') then
           if (i_pull > MXPULLING) then
              error_message = 'Error: too many PULL_UNRAVEL_CF in replica_pull (> MXPULLING), PROGRAM STOP'
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           i_pull = i_pull + 1
           !read(ctmp00, *) cdummy, imp, jmp, x, y, z
           read(ctmp00, *) cdummy, imp, jmp, inrep%pull_direction(1:3, i_pull)
           inmisc%ipull_unravel2mp(1,i_pull) = imp
           inmisc%ipull_unravel2mp(2,i_pull) = jmp
           !inmisc%pull_unravel_xyz(1,i_pull,:) = x
           !inmisc%pull_unravel_xyz(2,i_pull,:) = y
           !inmisc%pull_unravel_xyz(3,i_pull,:) = z
        endif
     end do ! iline
     
     ! checking input variables
     if (inrep%i_style(ivar) == REPVARSTYLE%LINEAR) then
        if (inrep%lowest(ivar) > INVALID_JUDGE) then
           error_message = 'Error: in reading value_lowest, PROGRAM STOP'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        if ( inrep%highest(ivar) > INVALID_JUDGE) then
           error_message = 'Error: in reading value_highest, PROGRAM STOP'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        
     else if (inrep%i_style(ivar) == REPVARSTYLE%EXPONENTIAL) then
        if ( inrep%lowest(ivar) > INVALID_JUDGE) then
           error_message = 'Error: in reading value_lowest, PROGRAM STOP'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        if ( inrep%highest(ivar) > INVALID_JUDGE) then
           error_message = 'Error: in reading value_highest, PROGRAM STOP'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        
     else if (inrep%i_style(ivar) == REPVARSTYLE%EXPLICIT) then
        do irep = 1, inrep%n_replica(ivar)
           if (inrep%var(irep,ivar) > INVALID_JUDGE) then
              write(error_message,*) 'Error: invalid value for REPLICA(', &
                   irep,') , PROGRAM STOP'
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        enddo
        
     else
        error_message = 'Error: invalid value for replica: i_style, PROGRAM STOP'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
     
     ! checking input variables
     if (ivar == REPTYPE%TEMP) then
        if (inrep%lowest(ivar) <= 0.0e0_PREC) then
           error_message = 'Error: invalid value for value_lowest of replica_temp, PROGRAM STOP'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
     endif
     if (ivar == REPTYPE%ION) then
        if (inrep%lowest(ivar) < 0.0e0_PREC) then
           error_message = 'Error: invalid value for value_lowest of replica_ion, PROGRAM STOP'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
     endif
     
  enddo ! do ivar = 1, REPTYPE%MAX

  inrep%n_pull = i_pull
  inmisc%npull_unravel = i_pull
  
#ifdef MPI_PAR
  end if
  call MPI_Bcast(n_replica_all,1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(n_dimension,  1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(flg_npar_rep, 1,          MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(flg_rep,      REPTYPE%MAX,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inrep,        inrep%sz,   MPI_BYTE,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inmisc, inmisc%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  
#endif
 

#ifdef _DEBUG
  write(6,*) ' setp_replica : n_replica_all = ' , n_replica_all
  if (inrep%flg_exchange) then
     write(6,*) ' setp_replica : i_exchange = 1'
  else
     write(6,*) ' setp_replica : i_exchange = 0'
  endif
  write(6,*) ' setp_replica : inrep%n_step_replica = ', inrep%n_step_replica
  write(6,*) ' setp_replica : inrep%n_step_save = '   , inrep%n_step_save
  write(6,*) ' setp_replica : inrep%npar_rep = '   , inrep%npar_rep
!  write(6,*) ' setp_replica : inrep%i_loadbalance = '   , inrep%i_loadbalance
!  write(6,*) ' setp_replica : inrep%n_step_adjust = '   , inrep%n_step_adjust
  do ivar = 1, REPTYPE%MAX
     write(6,*) ' setp_replica : DIMENSION',ivar
     write(6,*) ' setp_replica : inrep%i_style = '  , inrep%i_style(ivar)
     write(6,*) ' setp_replica : inrep%n_replica = '  , inrep%n_replica(ivar)
     write(6,*) ' setp_replica : inrep%lowest = '     , inrep%lowest(ivar)
     write(6,*) ' setp_replica : inrep%highest = '    , inrep%highest(ivar)
     !write(6,*) ' setp_replica : inrep%interval = '   , inrep%interval(ivar)
     do irep = 1, inrep%n_replica(ivar)
        write(6,*) ' setp_replica : inrep%var(',irep,') = ',inrep%var(irep,ivar)
     enddo
  enddo
  write(6,*) 'inp_replica_para : END'
#endif

endsubroutine inp_replica_para
