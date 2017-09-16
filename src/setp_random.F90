subroutine setp_random()

  use const_index
  use var_io, only: i_run_mode
  use var_setp, only: insimu, irand, mts, mts_rep
  use var_replica, only : n_replica_mpi, irep2grep
  use mt_stream
  use mt_kind_defs
  use mpiconst
#ifdef MPI_PAR
  use var_replica, only: n_replica_all
#endif

  implicit none
  integer       :: irep, grep, istream, tn
  integer(INT32) :: idstream, idstream_offset

  if(insimu%i_rand_type == 0) then
     allocate(mts(0:n_replica_mpi, 0:0))
  else
#ifdef MPI_PAR
     allocate(mts(0:n_replica_mpi*npar_mpi, 0:nthreads-1))
#else
     allocate(mts(0:n_replica_mpi, 0:nthreads-1))
#endif
  end if

  call set_mt19937
  call new(mts(0, 0))
  call init(mts(0, 0),irand)  ! init by array

  idstream_offset = 0

  if (i_run_mode == RUN%REPLICA) then 
     !! mts_rep is needed for consistency of random numbers between serial and parallel
     !! Because the internala state of mts(0,0) changes whenever create_stream is called,
     !! mts(0,0) may not be same between serial and parallel jobs.
     !! It causes different random number series in simu_replica_exchange(), 
     !! which is not desired particurarly for checking.
     call create_stream(mts(0, 0), mts_rep, 1)
     idstream_offset = 1
  endif     
  
  do irep = 1, n_replica_mpi

     if(insimu%i_rand_type == 0) then
        istream = irep
        grep = irep2grep(irep)
        idstream = grep
        if(idstream /= 0) then
           call create_stream(mts(0, 0), mts(istream, 0), idstream+idstream_offset)
        end if
     else
!$omp parallel private(tn, grep, istream, idstream)
        tn = 0
!$  tn = omp_get_thread_num()
        grep = irep2grep(irep)
        istream = irep
        idstream = (grep - 1)*nthreads + tn + 1
        if(idstream /= 0) then
           call create_stream(mts(0, 0), mts(istream, tn), idstream+idstream_offset)
        end if

#ifdef MPI_PAR
        if(insimu%i_rand_type == 1) then
           istream = irep + local_rank_mpi*n_replica_mpi
           idstream = ((grep - 1) + local_rank_mpi*n_replica_all)*nthreads + tn + 1
           if(idstream /= 0) then
              call create_stream(mts(0, 0), mts(istream, tn), idstream+idstream_offset)
           end if
        end if
#endif
!$omp end parallel

     end if
  end do
endsubroutine setp_random
