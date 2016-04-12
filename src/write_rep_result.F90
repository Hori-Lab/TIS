! write_rep_result
!> @brief Ouput for replica exchange simulation

subroutine write_rep_result()

  use const_maxsize
  use const_index
  use var_replica, only : n_replica_all, hist_attempt, hist_combination, &
                          hist_exchange, total_attempt, inrep, &
                          rate_exchange, flg_rep, n_turnover
  use var_inp,     only : outfile
#ifdef MPI_PAR
   use mpiconst
#endif

  implicit none

  integer :: rep_i, l_i, l_j, iperiod, ivar
  integer :: lunrep
  real(PREC) :: ratio(0:inrep%n_period_prob)
  integer :: total_hist_exchange
  integer :: total_hist_attempt
  integer :: n_sum_turnover(REPTYPE%MAX, n_replica_all)

#ifdef MPI_PAR
   if (myrank == 0) then
#endif

  lunrep = outfile%rep
  
  !#################################################
  ! Histogram
  !#################################################
  write(lunrep,*) ''
  write(lunrep,*) ''
  write(lunrep, '(a)') '## RESULT: Histogram of combination (LABEL:REPLICA)'
  write(lunrep, '(a6,1x)', ADVANCE='no') '#LABEL'
  do rep_i = 1, n_replica_all - 1
     write(lunrep, '(a3,i5.5,1x)', ADVANCE='no') 'REP', rep_i
  enddo
  write(lunrep, '(a3,i5.5)') 'REP', n_replica_all
  
  do l_i = 1, n_replica_all
     write(lunrep, '(i6,1x)',ADVANCE='no') l_i
     do rep_i = 1, n_replica_all - 1
        write(lunrep,'(i8,1x)', ADVANCE='no') hist_combination(l_i,rep_i)
     enddo
     write(lunrep,'(i8)') hist_combination(l_i,n_replica_all)
  enddo
  
  !#################################################
  ! Exchange rate
  !#################################################
  write(lunrep, *) ''
  write(lunrep, *) ''
  write(lunrep, '(a)') '## RESULT: Exchange rate'
  write(lunrep, '(a)', ADVANCE='no') '# Number of exchange attempt = '
  write(lunrep, *) total_attempt
  write(lunrep, '(a)') '# LABEL-LABEL  Rate'
  write(lunrep, '(a)', ADVANCE='no') '# Period          Total'
  do iperiod = 1, inrep%n_period_prob-1
     write(lunrep, '(1x,i8,2x)',ADVANCE='no') iperiod
  enddo
  write(lunrep, '(1x,i8)') inrep%n_period_prob

  do l_i = 1, n_replica_all -1
     loop_j: &
     do l_j = l_i+1, n_replica_all
        do iperiod = 1, inrep%n_period_prob
           if (hist_attempt(l_i, l_j, iperiod) == 0) then
              cycle loop_j
           endif
        enddo

        write(lunrep,'(2x,i5,a1,i5,1x)',ADVANCE='no') l_i,'-',l_j

        total_hist_exchange = 0
        total_hist_attempt = 0
        ratio(0:inrep%n_period_prob) = 0.0e0_PREC
        do iperiod = 1, inrep%n_period_prob
           total_hist_attempt  = total_hist_attempt  + hist_attempt(l_i, l_j, iperiod)
           ratio(0) = ratio(0) + real(hist_exchange(l_i, l_j, iperiod))
           ratio(iperiod) =  real(hist_exchange(l_i, l_j, iperiod)) &
                           / real(hist_attempt(l_i, l_j, iperiod))
        enddo

        ! total
        write(lunrep,'(1xf10.8)',ADVANCE='no')   ratio(0) / real(total_hist_attempt)
        ! period
        do iperiod = 1, inrep%n_period_prob-1
           write(lunrep,'(1xf10.8)', ADVANCE='no') ratio(iperiod)
        enddo
        write(lunrep,'(1xf10.8)') ratio(inrep%n_period_prob)
     enddo loop_j
  enddo

  !#################################################
  ! Average of exchange possibility
  !#################################################
  write(lunrep, *) ''
  write(lunrep, *) ''
  write(lunrep, '(a)') '## RESULT: Logarithmic average of exchange possibility' 
  write(lunrep, '(a)') '# LABEL-LABEL  Ave{ ln(Prob_exchange) }'
  write(lunrep, '(a)', ADVANCE='no') '# Period          Total'
  do iperiod = 1, inrep%n_period_prob-1
     write(lunrep, '(1x,i8,2x)',ADVANCE='no') iperiod
  enddo
  write(lunrep, '(1x,i8)') inrep%n_period_prob

  do l_i = 1, n_replica_all -1
     loop_j_2 : &
     do l_j = l_i+1, n_replica_all
        do iperiod = 1, inrep%n_period_prob
           if (hist_attempt(l_i, l_j, iperiod) == 0) then
              cycle loop_j_2
           endif
        enddo

        write(lunrep,'(2x,i5,a1,i5,1x)',ADVANCE='no') l_i,'-',l_j

        ratio(0:inrep%n_period_prob) = 0.0e0_PREC
        total_hist_attempt = 0
        do iperiod = 1, inrep%n_period_prob
           total_hist_attempt = total_hist_attempt + hist_attempt(l_i, l_j, iperiod)
           ratio(0) = ratio(0) + rate_exchange(l_i, l_j, iperiod)
           ratio(iperiod) =  rate_exchange(l_i, l_j, iperiod) &
                           / real(hist_attempt(l_i, l_j, iperiod), kind=PREC)
        enddo

        ! total
        write(lunrep,'(1xf10.4)',ADVANCE='no') ratio(0) / real(total_hist_attempt, kind=PREC)
        ! period
        do iperiod = 1, inrep%n_period_prob-1
           write(lunrep,'(1xf10.4)', ADVANCE='no') ratio(iperiod)
        enddo
        write(lunrep,'(1xf10.4)') ratio(inrep%n_period_prob)
     enddo loop_j_2
  enddo

  !#################################################
  ! Turnover
  !#################################################
  write(lunrep, *) ''
  write(lunrep, *) ''
  write(lunrep, '(a)') '## RESULT: Turnover number' 
  write(lunrep, '(a)', ADVANCE='no') '# Period                Total'
  do iperiod = 1, inrep%n_period_prob-1
     write(lunrep, '(1x,i14)',ADVANCE='no') iperiod
  enddo
  write(lunrep, '(1x,i14)') inrep%n_period_prob
  
  do ivar = 1, REPTYPE%MAX
     if (.not. flg_rep(ivar)) then
        cycle
     endif
     if (ivar == REPTYPE%TEMP) then
        write(lunrep, '(a)') '# Direction of temperature'
     elseif (ivar == REPTYPE%ION) then
        write(lunrep, '(a)') '# Direction of ionic strength'
     elseif (ivar == REPTYPE%PULL) then
        write(lunrep, '(a)') '# Direction of pull force'
     endif

     n_sum_turnover(:,:) = sum(n_turnover, 3)
     do rep_i = 1, n_replica_all
        write(lunrep, '(2x,a7,1x,i4)', ADVANCE='no') 'REPLICA', rep_i
        write(lunrep, '(1x, i14)', ADVANCE='no') n_sum_turnover(ivar, rep_i)
        do iperiod = 1, inrep%n_period_prob-1
           write(lunrep, '(1x,i14)',ADVANCE='no') n_turnover(ivar, rep_i, iperiod)
        enddo
        write(lunrep, '(1x,i14)')  n_turnover(ivar, rep_i, inrep%n_period_prob)
     enddo
  enddo

#ifdef MPI_PAR
   endif
#endif

endsubroutine write_rep_result
