! setp_anneal_para
!> @brief Read parameters for annealing simulation

! ***********************************************************************
subroutine setp_anneal_para()

  use const_maxsize
  use const_index
  use var_inp, only : infile, outfile
  use var_setp, only : inann

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  ! local variables
  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat

  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  ! --------------------------------------------------------------------
  ! default 
  inann%n_time_change  = -1
  inann%tempk_init = -1.0
  inann%tempk_last = -1.0

  ! --------------------------------------------------------------------

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'annealing       ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "annealing" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
            
        do iequa = 1, nequat
           cvalue = 'tempk_init'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inann%tempk_init, cvalue)
           
           cvalue = 'tempk_last'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inann%tempk_last, cvalue)
           
           cvalue = 'n_time_change'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inann%n_time_change, cvalue)
        end do
     end do
     
     ! -----------------------------------------------------------------
     ! checking input variables
     if(inann%tempk_init < 0.0) then
        error_message = 'Error: invalid value for tempk_init'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inann%tempk_last <= 0.0) then
        error_message = 'Error: invalid value for tempk_last'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inann%n_time_change <= 0.0) then
        error_message = 'Error: invalid value for n_time_change'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inann, inann%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_anneal_para
