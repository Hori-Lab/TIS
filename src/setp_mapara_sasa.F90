! setp_mapara_sasa

!> @brief Read parameters from sasa.para file. The parameters are &
!>        used for solvent accessible surface area dependent term

subroutine setp_mapara_sasa()
  use const_maxsize
  use const_index
  use var_inp, only: infile, outfile
  use var_setp, only : insasa

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! -------------------------------------------------------------------
  ! Local variables
  integer :: lunout, lunpara
  integer :: ioutput, iline, nlines
  integer :: aaid
  real(PREC) :: para, para1, para2, rsol, surf
  character(6) :: char6
  character(3) :: amino
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MSG_ERROR) :: error_message

  integer :: ifunc_seq2id

  ! -------------------------------------------------------------------
  ! Initialize parameters
  
  insasa%p_sasa(1:SEQID%MAX) = 0.0
  insasa%r_sasa(1:SEQID%MAX) = 0.0
  insasa%r_sol = 0.0
  insasa%coef_surf = 0.0

  ! -------------------------------------------------------------------
  ! Initialize file handle
  ioutput = 0
  lunpara = infile%para_fsasa
  lunout  = outfile%data

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  ! Reading input for P(i) and radius parameters
  rewind(lunpara)
     
  call ukoto_uiread2(lunpara, lunout, 'sasa_para        ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "sasa_para" in the sasa.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     read (ctmp00, *) amino, para1, para2
     aaid = ifunc_seq2id(amino)
     insasa%r_sasa(aaid) = para1
     insasa%p_sasa(aaid) = para2
  end do

  ! Reading input for connectivity
  rewind(lunpara)
  
  call ukoto_uiread2(lunpara, lunout, 'connectivity    ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "connectivity" in the sasa.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     read (ctmp00, *) char6, para
     if (char6 == 'RES1-2')insasa%connectivity(1) = para
     if (char6 == 'RES1-3')insasa%connectivity(2) = para
     if (char6 == 'RES1-4')insasa%connectivity(3) = para
     if (char6 == 'RES1-I')insasa%connectivity(4) = para
  end do
  
  ! Reading input for radius of solvent parameters
  rewind(lunpara)
  
  call ukoto_uiread2(lunpara, lunout, 'rsolvent        ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "rsolvent" in the sasa.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     read (ctmp00, *) rsol
     insasa%r_sol = rsol
  end do

  ! Reading input for coefficient of surface tension
  rewind(lunpara)

  call ukoto_uiread2(lunpara, lunout, 'surf_para       ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "surf_para" in the sasa.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     read (ctmp00, *) surf
     insasa%coef_surf = surf
  end do
#ifdef MPI_PAR
  end if
  call MPI_Bcast (insasa, insasa%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_mapara_sasa

