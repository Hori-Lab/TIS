! replica_runmode
!> @brief Read i_run_mode from .inp file

subroutine inp_runmode

  use const_maxsize
  use const_index
  use var_inp,     only : infile, i_run_mode

#ifdef MPI_PAR
  use mpiconst
#endif

  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat

  character(4)  :: kfind
  character(CARRAY_MXCOLM)  :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM)  :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  luninp = infile%inp
  lunout  = 6   ! At this moment, '.data' file is not exist.
  
  i_run_mode           = -1

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'job_cntl        ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)   
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "job_cntl" field in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     do iequa = 1, nequat
        cvalue = 'i_run_mode'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             i_run_mode, cvalue)
     end do
  end do
 
  if(i_run_mode == -1) then
     error_message = 'Error: cannot find "i_run_mode" in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  endif

#ifdef MPI_PAR
  end if
  call MPI_Bcast(i_run_mode,   1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine inp_runmode
