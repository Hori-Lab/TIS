!setp_cylinder_para
!> @brief Reads "cylinder_para" field and stores the data into  &
!>        "inmisc" struct.

subroutine setp_cylinder_para()

  use const_maxsize
  use const_index
  use const_physical
  use var_io, only : infile, outfile
  use var_setp, only : inmisc

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  ! intent(out) :: 

  ! --------------------------------------------------------------------
  ! local variables
  integer :: luninp, lunout
  integer :: iline, nlines, iequat, nequat

  character(4)  :: kfind
  character(CARRAY_MXCOLM)  :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MXCOLM)  :: cvalue
  character(CARRAY_MXCOLM)  :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data
  inmisc%cylinder_bgn  = 0.0
  inmisc%cylinder_end  = 0.0
  inmisc%cylinder_radi = 0.0
  inmisc%cylinder_coef = 0.0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'cylinder        ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "cylinder" field'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
          
     do iequat = 1, nequat
        cvalue = 'bgn'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             inmisc%cylinder_bgn, cvalue)
        cvalue = 'end'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             inmisc%cylinder_end, cvalue)
        cvalue = 'radi'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             inmisc%cylinder_radi, cvalue)
        cvalue = 'coef'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             inmisc%cylinder_coef, cvalue)
     enddo
  end do

  write(*, *) 'bgn: ', inmisc%cylinder_bgn
  write(*, *) 'end: ', inmisc%cylinder_end
  write(*, *) 'radi: ', inmisc%cylinder_radi
  write(*, *) 'coef: ', inmisc%cylinder_coef
  
  
#ifdef MPI_PAR
  end if

  call MPI_Bcast(inmisc, inmisc%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_cylinder_para
