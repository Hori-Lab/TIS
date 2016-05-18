! setp_para_ion
!> @brief Reads "initial_ion" block in input-file and stores the data into  &
!>        "inion" struct.

! ******************************************************************
subroutine setp_para_ion()
  
  use const_maxsize
  use const_index
  use var_io, only : infile, outfile
  use var_setp, only: inion

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  !------------------------------------------------------------------
  ! intent(out) :: 

  !------------------------------------------------------------------
  ! local variables
  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MXCOLM) :: ctmp02
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  ! ---------------------------------------------------------------------
  ! default setting 
  inion%num_na_ion = 0
  inion%num_k_ion = 0
  inion%num_cl_ion = 0
  inion%num_mg_ion = 0

  ! ---------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'initial_ion     ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "initial_ion" field in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
         
     do iequa = 1, nequat
        ctmp02 = csides(1, iequa)

        cvalue = 'num_na_ion'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inion%num_na_ion, cvalue)

        cvalue = 'num_k_ion'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inion%num_k_ion, cvalue)

        cvalue = 'num_cl_ion'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inion%num_cl_ion, cvalue)
        
        cvalue = 'num_mg_ion'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inion%num_mg_ion, cvalue)

     end do
  end do

  ! -----------------------------------------------------------------
  ! checking input variables
  if(inion%num_na_ion < 0) then
     error_message = 'Error: invalid value for num_na_ion'
     call util_error(ERROR%STOP_ALL, error_message)
  else
     write (lunout, *) "number of Na+ ion = ", inion%num_na_ion
  end if

  if(inion%num_k_ion < 0) then
     error_message = 'Error: invalid value for num_k_ion'
     call util_error(ERROR%STOP_ALL, error_message)
  else
     write (lunout, *) "number of K+ ion = ", inion%num_k_ion
  end if

  if(inion%num_cl_ion < 0) then
     error_message = 'Error: invalid value for num_cl_ion'
     call util_error(ERROR%STOP_ALL, error_message)
  else
     write (lunout, *) "number of Cl- ion = ", inion%num_cl_ion
  end if

  if(inion%num_mg_ion < 0) then
     error_message = 'Error: invalid value for num_mg_ion'
     call util_error(ERROR%STOP_ALL, error_message)
  else
     write (lunout, *) "number of Mg2+ ion = ", inion%num_mg_ion
  end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast (inion, inion%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

  inion%num_ion(IONTYPE%NA) = inion%num_na_ion
  inion%num_ion(IONTYPE%K)  = inion%num_k_ion
  inion%num_ion(IONTYPE%CL) = inion%num_cl_ion
  inion%num_ion(IONTYPE%MG) = inion%num_mg_ion

  inion%char_ion(IONTYPE%NA) = 'NA'
  inion%char_ion(IONTYPE%K ) = 'K '
  inion%char_ion(IONTYPE%CL) = 'CL'
  inion%char_ion(IONTYPE%MG) = 'MG'

end subroutine setp_para_ion
