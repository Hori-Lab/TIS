! setp_anchor_para
!> @brief This subroutine is to read and set the parameters for anchoring simulation, &
!>        in which some mass-points are constrained to some positions by harmonic springs.


! **********************************************************************
subroutine setp_anchor_para()

  use const_maxsize
  use const_index
  use var_io, only : infile, outfile
  use var_setp, only : inmisc

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  ! intent(out) :: nanc, ianc2mp, coef_anc, anc_dist, anc_xyz

  ! --------------------------------------------------------------------
  ! local variables
  integer :: luninp, lunout
  integer :: iline, nlines
  integer :: ianc, ianc_com_ini
  character(4) :: char4
  character(15) :: char15
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data
  inmisc%nanc = 0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'anchor_para     ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "anchor_para" field in set_anchor_para'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  ianc = 0
  ianc_com_ini = 0
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     if (ctmp00(1:15) == 'ANCH_CENTER_INI') then
        ianc_com_ini = ianc_com_ini + 1
        read (ctmp00, *) &
             char15, &
             inmisc%ianc_com_ini2grp(ianc_com_ini), &
             inmisc%coef_anc_com_ini(ianc_com_ini), &
             inmisc%anc_com_ini_dist(ianc_com_ini)
        write (lunout, '(2a, i10)') &
             '---reading anchor group (initial position): ', &
             char15, &
             inmisc%ianc_com_ini2grp(ianc_com_ini), &
             inmisc%coef_anc_com_ini(ianc_com_ini), &
             inmisc%anc_com_ini_dist(ianc_com_ini)
     else if(ctmp00(1:4) == 'ANCH') then
        ianc = ianc + 1
        read (ctmp00, *) char4, inmisc%ianc2mp(ianc), &
             inmisc%coef_anc(ianc), inmisc%anc_dist(ianc), &
             inmisc%anc_xyz(1, ianc), inmisc%anc_xyz(2, ianc), &
             inmisc%anc_xyz(3, ianc)
        write (lunout, '(2a, i10, 5g10.3)') '---reading anchor residue: ', &
             char4, inmisc%ianc2mp(ianc), &
             inmisc%coef_anc(ianc), inmisc%anc_dist(ianc), &
             inmisc%anc_xyz(1, ianc), inmisc%anc_xyz(2, ianc), &
             inmisc%anc_xyz(3, ianc)
     end if
  end do
  inmisc%nanc = ianc
  inmisc%nanc_com_ini = ianc_com_ini
  
  if(inmisc%nanc > MXANCHOR) then
     error_message = 'Error: should increase MXANCHOR'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inmisc, inmisc%sz ,MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_anchor_para
