! setp_box
!> @brief This subroutine is to read and set the parameters for performing the simulation with a rectangular box.

subroutine setp_box()

  use const_maxsize
  use const_index
  use var_inp, only : infile, outfile
  use var_setp, only : inmisc

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  ! intent(out) :: xbox, ybox, zbox

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
  inmisc%xbox = -1.0
  inmisc%ybox = -1.0
  inmisc%zbox = -1.0
  inmisc%boxsigma = -1.0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'in_box          ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "in_box" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

   
     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
   
        do iequa = 1, nequat
           cvalue = 'xbox'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmisc%xbox, cvalue)

           cvalue = 'ybox'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmisc%ybox, cvalue)

           cvalue = 'zbox'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmisc%zbox, cvalue)

           cvalue = 'boxsigma'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmisc%boxsigma, cvalue)
        end do
     end do

     ! -----------------------------------------------------------------
     ! checking input variables
     if(inmisc%xbox <= 0.0) then
        error_message = 'Error: invalid value for xbox'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inmisc%ybox <= 0.0) then
        error_message = 'Error: invalid value for ybox'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inmisc%zbox <= 0.0) then
        error_message = 'Error: invalid value for zbox'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inmisc%boxsigma <= 0.0) then
        error_message = 'Error: invalid value for boxsigma'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inmisc, inmisc%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)

#endif

end subroutine setp_box


