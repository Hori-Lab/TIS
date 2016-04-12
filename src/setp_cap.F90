! setp_box
!> @brief This subroutine is to read and set the parameters for performing the simulation with a spherical cap

subroutine setp_cap()

  use const_maxsize
  use const_index
  use var_inp, only : infile, outfile
  use var_setp, only : inmisc

#ifdef MPI_PAR
  use mpiconst
#endif
  
  implicit none

  ! --------------------------------------------------------------------
  ! local variables
  integer :: i ! tmp pointer for do loop
  integer :: luninp, lunout
  integer :: nlines, iline, iequa, nequat
  
  character(4) :: kfind
  character(6) :: char6
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue, ctmp00
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message
  
  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data
  inmisc%rcap = 1.0e2_PREC ! default value is '100.0'
  inmisc%kcap = 5.0e2_PREC ! default value is '500.0'
  do i = 1, 3
     ! default value is (0.0, 0.0, 0.0)
     inmisc%center_cap(i) = 0.0e0_PREC 
  end do

#ifdef MPI_PAR
  if (myrank == 0) then
#endif
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'in_cap          ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "in_cap" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        ctmp00 = cwkinp(iline)
        
        call ukoto_uiequa2(lunout, ctmp00, nequat, csides)

        do iequa = 1, nequat
           cvalue = 'rcap'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmisc%rcap, cvalue)
           cvalue = 'kcap'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmisc%kcap, cvalue)
        end do

        if (ctmp00(1:6) == 'CENTER') then
           read(ctmp00, *) char6, inmisc%center_cap(1), inmisc%center_cap(2), inmisc%center_cap(3)
        end if
        
     end do

#ifdef MPI_PAR
  end if
  
  call MPI_Bcast(inmisc, inmisc%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif
  
  
end subroutine setp_cap
