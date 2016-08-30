subroutine setp_random()

  use var_setp, only: insimu, irand, mts
  use var_replica, only : n_replica_mpi, irep2grep
  use mt_stream
  use mt_kind_defs
  use mpiconst
#ifdef MPI_PAR
  use var_replica, only: n_replica_all
#endif

  implicit none
  integer       :: irep, grep, istream, tn
  integer(INT32) :: idstream

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
  
  do irep = 1, n_replica_mpi

     if(insimu%i_rand_type == 0) then
        istream = irep
        grep = irep2grep(irep)
        idstream = grep
        if(idstream /= 0) then
           call create_stream(mts(0, 0), mts(istream, 0), idstream)
        end if
     else
!$omp parallel private(tn, grep, istream, idstream)
        tn = 0
!$  tn = omp_get_thread_num()
        grep = irep2grep(irep)
        istream = irep
        idstream = (grep - 1)*nthreads + tn + 1
        if(idstream /= 0) then
           call create_stream(mts(0, 0), mts(istream, tn), idstream)
        end if

#ifdef MPI_PAR
        if(insimu%i_rand_type == 1) then
           istream = irep + local_rank_mpi*n_replica_mpi
           idstream = ((grep - 1) + local_rank_mpi*n_replica_all)*nthreads + tn + 1
           if(idstream /= 0) then
              call create_stream(mts(0, 0), mts(istream, tn), idstream)
           end if
        end if
#endif
!$omp end parallel

     end if
  end do
endsubroutine setp_random
