subroutine setp_twobody_gauss()

  use const_maxsize
  use const_index
  use const_physical
  use var_io, only : infile, outfile
  use var_setp, only : inmisc
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

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

  inmisc%exv_gauss_a0 = INVALID_VALUE
  inmisc%con_gauss_sigma = INVALID_VALUE
  inmisc%con_gauss_k = INVALID_VALUE

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'twobody_gauss   ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "twobody_gauss" field in set_twobody_gauss'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)

     do iequa = 1, nequat
        cvalue = 'exv_gauss_a0'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inmisc%exv_gauss_a0, cvalue)

        cvalue = 'con_gauss_sigma'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inmisc%con_gauss_sigma, cvalue)

        cvalue = 'con_gauss_k'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inmisc%con_gauss_k, cvalue)

     end do
  end do

  ! -----------------------------------------------------------------
  ! checking input variables
  if(inmisc%exv_gauss_a0 > INVALID_JUDGE) then
     error_message = 'Error: invalid value for exv_gauss_a0 in twobody_gauss field'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif(inmisc%con_gauss_sigma > INVALID_JUDGE) then
     error_message = 'Error: invalid value for con_gauss_sigma in twobody_gauss field'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif(inmisc%con_gauss_k > INVALID_JUDGE) then
     error_message = 'Error: invalid value for con_gauss_k in twobody_gauss field'
     call util_error(ERROR%STOP_ALL, error_message)

  endif

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inmisc, inmisc%sz ,MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_twobody_gauss
