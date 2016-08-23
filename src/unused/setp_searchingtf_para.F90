! setp_searchingtf_para
!> @brief This subroutine is to read and set the parameters for the automatic T_F (folding temperature) search.

subroutine setp_searchingtf_para()
  
  use const_maxsize
  use const_index
  use var_io, only : infile, outfile
  use var_setp, only : insear
  use var_struct, only : nunit_real, nmp_real, iclass_unit

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  insear%tempk_upper = -1.0
  insear%tempk_lower = -1.0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif
   
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'searching_tf    ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
   
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "searching_tf" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
        do iequa = 1, nequat

           cvalue = 'tempk_upper'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                insear%tempk_upper, cvalue)

           cvalue = 'tempk_lower'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                insear%tempk_lower, cvalue)

        end do
     end do
   
     ! -----------------------------------------------------------------
     ! checking input variables
     if(insear%tempk_upper < 0.0) then
        error_message = 'Error: invalid value for tempk_upper'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(insear%tempk_lower < 0.0) then
        error_message = 'Error: invalid value for tempk_lower'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     if(nunit_real > 1) then
        error_message = 'Error: Seaching Tf is not available for more than 1 chain'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(nmp_real > 150) then
        error_message = 'Error: Seaching Tf is not available for more than 150 amino acid'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(iclass_unit(1) /= CLASS%PRO) then
        error_message = 'Error: Seaching Tf is available only for protein'
        call util_error(ERROR%STOP_ALL, error_message)

     end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast(insear, insear%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif
  
end subroutine setp_searchingtf_para
