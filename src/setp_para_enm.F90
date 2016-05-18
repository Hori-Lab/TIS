!setp_para_enm
!> @brief Reads "elastic_network" field and stores the data into  &
!>        "inenm" struct.

subroutine setp_para_enm()

  use const_maxsize
  use const_index
  use var_io,   only : infile, outfile
  use var_setp, only : inenm
  use mpiconst

  implicit none

  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat

  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! ----------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'elastic_network  ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "elastic_network" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
   
     ! ----------------------------------------------------------------------
     ! reading parameters
     ! default parameters
     inenm%cenm = -1.0
     inenm%dfcontact_enm = -1.0
     
     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
            
        do iequa = 1, nequat
           cvalue = 'cenm'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inenm%cenm, cvalue)
           
           cvalue = 'dfcontact_enm'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inenm%dfcontact_enm, cvalue)
        end do
     end do

     ! -----------------------------------------------------------------
     ! checking input variables
     if(inenm%cenm <= 0.0) then
        error_message = 'Error: invalid value for cenm'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inenm%dfcontact_enm <= 0.0) then
        error_message = 'Error: invalid value for dfcontact_enm'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inenm%cenm,          1,PREC_MPI   ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inenm%dfcontact_enm, 1,PREC_MPI   ,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_para_enm
