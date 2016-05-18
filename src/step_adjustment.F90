#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

!#define _DEBUG_ADJUST

subroutine step_adjustment(istep, n_exchange, i_loadbalance)
  use const_maxsize
  use const_index
  use var_setp,    only : inmisc, insimu
  use var_struct,  only : lele
  use var_replica, only : inrep, n_replica_all, n_replica_mpi, &
                          irep2grep, lab2step, rep2step, exchange_step
  use var_io,     only : outfile
  use time
#ifdef MPI_PAR
  use mpiconst
#endif
  implicit none
  ! ------------------------------------------------
  integer(L_INT), intent(in) :: istep 
  integer, intent(in) :: n_exchange
  integer, intent(in) :: i_loadbalance
  ! ------------------------------------------------
  integer    :: irep, grep
  integer    :: nstep, maxstep
  integer    :: lele_max
  integer    :: lele_inn(MXREPLICA)
  integer    :: lele_all(MXREPLICA)  ! number of neighbor list for electrostatic interaction
  integer    :: lunout
  real(PREC) :: loop_time, loop_time_rep(MXREPLICA)
  real(PREC) :: loop_time_max_rep
  !real(PREC) :: loop_time_min_rep
  real(PREC) :: correct_step, ratio
  parameter (ratio = 1.00_PREC)
  integer, save :: n_adjust = 0
#ifdef MPI_PAR
  real(PREC) :: loop_time_max_inn
#endif

  lunout = outfile%data
!
! maximum number of exchange step
!
  nstep   = inrep%n_step_replica
  maxstep = 10 * nstep

  write(6,'(a,i3,i10)') ' <<< step_adjustment was called (i_loadbalance,istep)>>>', i_loadbalance,istep

  if (i_loadbalance == 1) then
!
! balanced by neighborlist_ele
!
     if (inmisc%force_flag(INTERACT%ELE)) then
        lele_inn(:) = 0
        lele_all(:) = 0

        TIME_S(tmc_step_adj)
#ifdef MPI_PAR

#ifdef SHARE_NEIGH
!
! comminication for replica parallel process
!
        call mpi_allgather(lele    , n_replica_mpi, MPI_INTEGER, &
                           lele_all, n_replica_mpi, MPI_INTEGER, &
                           mpi_comm_rep,ierr)
#else
!
! comminication for inner parallel process
!
        call mpi_allreduce(lele, lele_inn, n_replica_mpi, MPI_INTEGER, &
                           MPI_SUM, mpi_comm_local, ierr)

!       call mpi_allreduce(lele, lele_inn, n_replica_mpi, MPI_INTEGER, &
!                          MPI_SUM, mpi_comm_local, ierr)

!
! comminication for replica parallel process
!
        call mpi_allgather(lele_inn, n_replica_mpi, MPI_INTEGER, &
                           lele_all, n_replica_mpi, MPI_INTEGER, &
                           mpi_comm_rep,ierr)
#endif

#else
        lele_all(1:n_replica_mpi) = lele(1:n_replica_mpi)
#endif
        TIME_E(tmc_step_adj)

        lele_max = maxval(lele_all(:)) 

        do irep = 1, n_replica_all
           if (inrep%n_step_replica > 0) then
!            write(6,*) ' = irep, lele_max,lele_all(irep) ', irep, lele_max, lele_all(irep)
             if (lele_all(irep) /= 0) then
                correct_step = dble((lele_max-lele_all(irep))) / dble((lele_all(irep))) * ratio
                lab2step(irep) = min(maxstep, int(nstep*(1.0_PREC+correct_step)))
!               lab2step(irep) = min(maxstep, (lele_max)*nstep/(lele_all(irep)))
!               lab2step(irep) = nstep + (irep-1)*nstep/n_replica_all
             endif
           else
              lab2step(irep) = 1
           endif
        enddo

#ifdef _DEBUG_ADJUST
! number of neighborlist for electrostatic interaction
        do irep=1, n_replica_all
           write(6,'(a,2i6)') 'step_adjustment,irep,lele_all(irep) ' ,irep,lele_all(irep)
        enddo
        call flush(6)
#endif

     endif  ! if (inmisc%force_flag(INTERACT%ELE))
!
  else if (i_loadbalance == 2) then
!
! balanced by measurment time
!
     if (istep == insimu%i_tstep_init) then
        return
     endif
     loop_time = total_time(tm_lap)
#ifdef _DEBUG_ADJUST
     write(6,'(a18,f15.6)') ' main_loop time' ,loop_time
#endif

     TIME_S(tmc_step_adj)
#ifdef MPI_PAR
     call mpi_allreduce(loop_time, loop_time_max_inn, 1, PREC_MPI, &
                        MPI_MAX, mpi_comm_local, ierr)

     call mpi_allgather(loop_time_max_inn, 1, PREC_MPI, &
                        loop_time_rep,     1, PREC_MPI, &
                        mpi_comm_rep, ierr)
!
!    call mpi_allreduce(loop_time, loop_time_max_inn, 1, PREC_MPI, &
!                       MPI_MAX, mpi_comm_local, ierr)
!    call mpi_allreduce(loop_time_max_inn, loop_time_max_rep, 1, PREC_MPI, &
!                       MPI_MAX, mpi_comm_rep, ierr)
!    call mpi_allreduce(loop_time_max_inn, loop_time_min_rep, 1, PREC_MPI, &
!                       MPI_MIN, mpi_comm_rep, ierr)
#ifdef _DEBUG_ADJUST
     write(6,'(a22,f15.6)') 'loop_time_max_inn ' ,loop_time_max_inn
#endif
#endif
     TIME_E(tmc_step_adj)

#ifdef _DEBUG_ADJUST
     do irep=1, n_replica_all
        write(6,'(a,i5,f15.6)') 'irep, loop_time_rep ' ,irep, loop_time_rep(irep)
     enddo
#endif
!
     loop_time_max_rep = maxval(loop_time_rep(1:n_replica_all)) 
     do irep = 1, n_replica_all
        if (inrep%n_step_replica > 0) then
          correct_step = (loop_time_max_rep-loop_time_rep(irep)) / (loop_time_rep(irep)) * ratio
          lab2step(irep) = min(maxstep, int(nstep*(1.0_PREC+correct_step)))

!         lab2step(irep) = min(maxstep, int((loop_time_max_rep)*nstep/(loop_time_rep(irep))))
!         lab2step(irep) = min(maxstep, int((loop_time_max_rep)*nstep/(loop_time)))
        else
          lab2step(irep) = 1
        endif
     enddo

  endif  ! ------  if (i_loadbalance)

  if (istep == insimu%i_tstep_init) then
     do irep=1, n_replica_mpi
        grep = irep2grep(irep)
        exchange_step(grep) = rep2step(grep)
#ifdef _DEBUG_ADJUST
        write(6,'(a,4i5)') ' - irep,grep,exchange_step(grep),rep2step(grep) - ' ,\
                               irep,grep,exchange_step(grep),rep2step(grep)
#endif
     enddo
  endif
!
#ifdef _DEBUG_ADJUST
  do irep = 1, n_replica_all
     write(6,*) ' = irep, lab2step(irep) = ', irep, lab2step(irep)
  enddo
#endif

  n_adjust = n_adjust + 1
  call write_rep_adjust(istep, n_exchange, n_adjust, i_loadbalance, lele_all, loop_time_rep)

  return
end subroutine step_adjustment
