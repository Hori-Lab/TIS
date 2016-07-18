! setp_fix_para
!> @brief This subroutine is to read and set the parameters &
!>        to fix some mass-points or some units to their initial positions.

! **********************************************************************
subroutine setp_fix_para()

  use const_maxsize
  use const_index
  use var_io, only : infile, outfile
  use var_setp, only : fix_mp
  use var_struct, only : lunit2mp

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  ! local variables
  integer :: imp, iunit, icol
  integer :: luninp, lunout
  integer :: iline, nlines
  integer :: inunit(2), instate

  character(12) :: char12
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  fix_mp(:) = .False.

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'fix_para        ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "fix_para" filed in set_fix_para'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
   
     if(ctmp00(1:8) == 'FIX_UNIT') then
        do icol = 10, CARRAY_MXCOLM
           if(ctmp00(icol:icol) == ')') exit
        end do
        
        read (ctmp00(10:icol-1), *) char12
        write (lunout, '(2a)') '---reading FIX_UNIT: ', trim(char12)
        call util_unitstate(char12, inunit, instate)

        do iunit = inunit(1), inunit(2)
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
              fix_mp(imp) = .True.
           end do
        end do

     else if(ctmp00(1:6) == 'FIX_MP') then
        do icol = 8, CARRAY_MXCOLM
           if(ctmp00(icol:icol) == ')') exit
        end do
        
        read (ctmp00(8:icol-1), *) char12
        write (lunout, '(2a)') '---reading FIX_MP: ', trim(char12)
        call util_unitstate(char12, inunit, instate)

        do imp = inunit(1), inunit(2)
           fix_mp(imp) = .True.
        end do

     end if
     
  end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast(fix_mp, nmp_all, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

end subroutine setp_fix_para
