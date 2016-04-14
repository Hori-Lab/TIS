!input
!> @brief Receives the name of inputfile given to CafeMol, and calls  &
!>        "inp_replica", "inp_filename", "inp_job" "inp_datafile".

!#define _DEBUG
subroutine input()

  use const_maxsize
  use const_index
  use var_inp,     only : outfile, infile, i_run_mode, flg_rst
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  !----------------------------------------------------------------------
  ! local variables
  integer :: n
  integer :: iopen_status
  character(CARRAY_MSG_ERROR) :: error_message
#ifdef _DEBUG
  integer :: i_debug
#endif

  ! For getting inputfile's name
  integer :: iarg                  ! # of arguments excepting the program name
  integer :: iargc
  character(CARRAY_MXFILE) :: cfile_input    ! input filename

  !----------------------------------------------------------------------
  iopen_status = 1
  iarg         = 0
  flg_rst = .false.

  !----------------------------------------------------------------------
  ! read input file and store it in scratch file

#ifdef _DEBUG
  write(*,*) '#### start input'
#endif

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  iarg = iargc()
  ! exception
  if (iarg < 1 .OR. iarg > 2) then
     error_message = 'Usage: % PROGRAM (INPUT_FILE) [(restart file)]'
     call util_error(ERROR%STOP_STD, error_message)
  end if

  !==============
  ! input file
  !==============
  call getarg(1, cfile_input) 
  open (infile%inp, file=cfile_input, status='OLD', action='READ', iostat = iopen_status)
  ! exception
  if(iopen_status > 0) then
     error_message = 'Error: cannot open the file in input '//cfile_input
     call util_error(ERROR%STOP_STD, error_message)
  end if

  !==============
  ! restart file
  !==============
  if (iarg == 2) then
     call getarg(2, cfile_input)
     open(infile%rst, file=cfile_input, status='old', action='read', &
          iostat=iopen_status, &
#ifdef UNFORMATTED
          form='unformatted', access='transparent') 
#else
!          form='binary')
          form = 'unformatted', access = 'stream')
#endif
     ! exception
     if(iopen_status > 0) then
        error_message = 'Error: cannot open the file in input '//cfile_input
        call util_error(ERROR%STOP_STD, error_message)
     end if
     n = index(cfile_input,' ')
     write (*, *) "open restart file:", cfile_input(1:n-1)
     flg_rst = .true.
  endif

#ifdef MPI_PAR
  end if
  call MPI_Bcast(flg_rst,1,MPI_LOGICAL,0 ,MPI_COMM_WORLD,ierr)
#endif

  !----------------------------------------------------------------------
  ! reading i_run_mode for replica exchange
  call inp_runmode()

  write (*, *) "i_run_mode = ", i_run_mode

  ! reading parameters for replica exchange
  call inp_replica_para()


  ! parallelization settings for MPI 
  call inp_split_mpi()

  
  ! making the replica set tables
  if (i_run_mode == RUN%REPLICA) then
     call inp_replica_tables()
  end if


  ! open data file and movie file
  call inp_filename()

  ! reading parameters of job_control
  call inp_job()

  ! reading unit_and_state
  call inp_unitstate()

  ! -----------------------------------------------------------------------
  ! reading energy function
  call inp_energy_func()
  
  ! open parameter files
  call inp_datafile()

  ! output replica variable table to .rep file
  if (i_run_mode == RUN%REPLICA) then
     call write_rep_table()
  endif

#ifdef _DEBUG
  write(6,*) 'input : infile%inp = ',infile%inp
  write(6,*) 'input : outfile%data = ',outfile%data
  write(6,*) 'input : outfile%ninfo = ',outfile%ninfo
  do i_debug = 1, MXREPLICA
     write(6,*) 'input : outfile%ts(',i_debug,') = ',outfile%ts(i_debug)
  enddo
  do i_debug = 1, MXREPLICA
     write(6,*) 'input : outfile%movie(',i_debug,') = ',outfile%movie(i_debug)
  enddo
  do i_debug = 1, MXREPLICA
     write(6,*) 'input : outfile%crd(',i_debug,') = ',outfile%crd(i_debug)
  enddo
  do i_debug = 1, MXREPLICA
     write(6,*) 'input : outfile%velo(',i_debug,') = ',outfile%velo(i_debug)
  enddo
  do i_debug = 1, MXREPLICA
     write(6,*) 'input : outfile%dcd(',i_debug,') = ',outfile%dcd(i_debug)
  enddo
  do i_debug = 1, MXREPLICA
     write(6,*) 'input : outfile%vdcd(',i_debug,') = ',outfile%vdcd(i_debug)
  enddo
  do i_debug = 1, MXREPLICA
     write(6,*) 'input : outfile%pdb(',i_debug,') = ',outfile%pdb(i_debug)
  enddo
  write(6,*) 'input : outfile%rep = ',outfile%rep
  write(*,*) '#### end input'
#endif

end subroutine input
