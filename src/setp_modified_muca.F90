!setp_modified_muca
!> @brief Reads "<<<< modified_muca" field. Data are stored in  &
!>        "inmmc" struct.

subroutine setp_modified_muca()

  use const_maxsize
  use const_index
  use var_io, only : infile, outfile
  use var_setp, only : inmmc

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
  inmmc%em_depth = -1.0
  inmmc%em_mid = -1.0
  inmmc%em_sigma = -1.0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'modified_muca     ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "modified_muca" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

   
     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
   
        do iequa = 1, nequat
           cvalue = 'em_depth'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmmc%em_depth, cvalue)

           cvalue = 'em_mid'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmmc%em_mid, cvalue)

           cvalue = 'em_sigma'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmmc%em_sigma, cvalue)
        end do
     end do

     ! -----------------------------------------------------------------
     ! checking input variables
     if(inmmc%em_depth <= 0.0) then
        error_message = 'Error: invalid value for em_depth'
        call util_error(ERROR%STOP_ALL, error_message)

!     else if(inmmc%em_mid <= 0.0) then
!        error_message = 'Error: invalid value for em_mid'
!        call util_error(ERROR%STOP_ALL, error_message)

     else if(inmmc%em_sigma <= 0.0) then
        error_message = 'Error: invalid value for em_sigma'
        call util_error(ERROR%STOP_ALL, error_message)

     end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inmmc, inmmc%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)

#endif

end subroutine setp_modified_muca


