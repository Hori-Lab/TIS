! inp_split_mpi
!> @brief Constructs the parallelization settings for MPI  &
!>        in REMD simulation.

subroutine inp_split_mpi()

  use const_maxsize
  use const_index
  use var_io,     only : i_run_mode
  use var_replica, only : inrep, n_replica_all, flg_npar_rep, &
                          n_replica_mpi, irep2grep, grep2irep, grep2rank
  use mpiconst

  implicit none

  integer :: irep
#ifdef MPI_PAR
  integer :: iproc, grep
  integer :: jsta, jend, jlen
  character(CARRAY_MSG_ERROR) :: error_message
#endif

#ifdef MPI_PAR

  if(i_run_mode /= RUN%REPLICA .OR. nprocs == 1) then
     npar_rep = 1
  else
     if(flg_npar_rep) then
        npar_rep = inrep%npar_rep
     else
        if(nprocs >= n_replica_all) then
           npar_rep = n_replica_all
        else
           npar_rep = nprocs
        end if
        inrep%npar_rep = npar_rep
     end if
  end if

  if (nprocs ==  1) then
    npar_tmp = 1
  else
    npar_tmp = nprocs / npar_rep   ! npar_tmp = npar_mpi
    if (nprocs .ne. npar_tmp * npar_rep) then
      write(error_message,*) "Error: nprocs can't devide by npar_rep "
      call util_error(ERROR%STOP_ALL, error_message)
!      write(6,*) " warning nprocs can't devide by npar_rep "
!      do while (mod(nprocs, npar_rep) .ne. 0)
!        npar_rep = npar_rep -1
!        write(6,*) ' npar_rep ', npar_rep
!      end do
!      npar_tmp = nprocs / npar_rep
!      write(6,*) " procs, npar_rep, npar_mpi ", nprocs, npar_rep, npar_tmp
    end if
  end if

  icolor = myrank / npar_tmp
  ikey   = myrank

  call MPI_barrier(MPI_COMM_WORLD,ierr)

  call flush(6)

  call mpi_comm_split(MPI_COMM_WORLD, icolor, ikey, mpi_comm_local,ierr)

  call mpi_comm_size(mpi_comm_local, npar_mpi, ierr)
  call mpi_comm_rank(mpi_comm_local, local_rank_mpi, ierr)

  icolor = mod(myrank, npar_mpi)
  ikey   = myrank
 
  call mpi_comm_split(MPI_COMM_WORLD, icolor, ikey, mpi_comm_rep,ierr)
  call mpi_comm_size(mpi_comm_rep, npar_rep, ierr)
  call mpi_comm_rank(mpi_comm_rep, local_rank_rep, ierr)

  grep2rank(:) = 0
  grep2irep(:) = 0
  do iproc = 0, nprocs-1
     if (myrank == iproc) then

        print *,"[mpi] npar_mpi       = ",npar_mpi
        print *,"[mpi] local_rank_mpi = ",local_rank_mpi
      
        print *,"[mpi] npar_rep       = ",npar_rep
        print *,"[mpi] local_rank_rep = ",local_rank_rep
      
        jlen = (n_replica_all-1+npar_rep)/npar_rep
        jsta = 1+jlen*local_rank_rep
        jend = min(jsta+jlen-1, n_replica_all)
        n_replica_mpi = jend - jsta + 1
        print *,"[mpi] n_replica_mpi  = ",n_replica_mpi
        call flush(6)
      
        do irep = 1, n_replica_mpi
           grep = jlen*local_rank_rep + irep
           irep2grep(irep) = grep
           grep2irep(grep) = irep
           grep2rank(grep) = local_rank_rep
           print *,"[mpi] irep2grep(",irep,") = ",irep2grep(irep)
        enddo

! write(6,'(a,5i5)') ' myrank, icolor, ikey, nprocs_local, local_rank ', &
!                      myrank, icolor, ikey, nprocs_local, local_rank
        call flush(6)
   
     endif
     call MPI_barrier(MPI_COMM_WORLD,ierr)
  enddo

  call MPI_allreduce(MPI_IN_PLACE, grep2rank, MXREPLICA, MPI_INTEGER, &
                     MPI_SUM, mpi_comm_rep, ierr)
  call MPI_allreduce(MPI_IN_PLACE, grep2irep, MXREPLICA, MPI_INTEGER, &
                     MPI_SUM, mpi_comm_rep, ierr)
  do irep = 1, n_replica_all
     print *,"[mpi] grep2irep(",irep,") = ",grep2irep(irep)
     print *,"[mpi] grep2rank(",irep,") = ",grep2rank(irep)
  enddo

! ---------------------------------------------------------------------
! Not MPI
#else
  local_rank_mpi = 0
  local_rank_rep = 0
  n_replica_mpi = n_replica_all
  print *,"[single] n_replica_mpi  = ",n_replica_mpi
  do irep = 1, n_replica_all
     irep2grep(irep) = irep
     grep2irep(irep) = irep
     grep2rank(irep) = 0
     print *,"[single] irep2grep(",irep,") = ",irep2grep(irep)
     print *,"[single] grep2irep(",irep,") = ",irep2grep(irep)
  enddo

#endif

endsubroutine inp_split_mpi
