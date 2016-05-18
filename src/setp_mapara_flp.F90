! setp_mapara_lfp

!> @brief Read parameters from flexible_local.para file. The parameters are &
!>        used for flexible local interaction

subroutine setp_mapara_flp()
  use const_maxsize
  use const_index
  use var_io, only: infile, outfile
  use var_setp, only : inflp

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! -------------------------------------------------------------------
  ! Local variables
  integer :: lunout, lunpara
  integer :: ioutput, iline, nlines
  integer :: ipara, npara_dih, npara_ang
  integer :: aaid(2)
  real(PREC) :: param(10), x

  character(3) :: amino1, amino2
  character(4) :: para
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MSG_ERROR) :: error_message

  integer :: ifunc_seq2id

  ! -------------------------------------------------------------------
  ! Initialize parameters
  
  inflp%ang_para_x(1:10) = -1.0
  inflp%ang_para_y(1:20, 1:10) = -1.0
  inflp%ang_para_y2(1:20, 1:10) = -1.0
  inflp%dih_para(1:400, 1:400, 1:7) = -1.0
  
  ! -------------------------------------------------------------------
  ! Initialize file handle
  ioutput = 0
  lunpara = infile%para_flp
  lunout  = outfile%data
  npara_dih = 7
  npara_ang = 10
  
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  ! Reading input for dihedral angle parameters
  rewind(lunpara)
     
  call ukoto_uiread2(lunpara, lunout, 'dihedral_angle  ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "dihedral_angle" in the flexible_local.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     read (ctmp00, *) amino1, amino2, & 
          param(1), param(2), param(3), param(4), param(5), param(6), param(7)
     aaid(1) = ifunc_seq2id(amino1)
     aaid(2) = ifunc_seq2id(amino2)
     do ipara = 1, npara_dih
        inflp%dih_para(aaid(1), aaid(2), ipara) = param(ipara)
     end do
  end do

  ! Reading input for bond angle parameters about x
  rewind(lunpara)
  
  call ukoto_uiread2(lunpara, lunout, 'bond_angle_x    ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "bond_angle_x" in the flexible_local.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     read (ctmp00, *) para, x
     inflp%ang_para_x(iline) = x
  end do
  
  ! Reading input for bond angle parameters about y
  rewind(lunpara)
  
  call ukoto_uiread2(lunpara, lunout, 'bond_angle_y    ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "bond_angle_y" in the flexible_local.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     read (ctmp00, *) amino1, &
          param(1), param(2), param(3), param(4), param(5), &
          &param(6), param(7), param(8), param(9), param(10)
     aaid(1) = ifunc_seq2id(amino1)
     
     do ipara = 1, npara_ang
        inflp%ang_para_y(aaid(1), ipara) = param(ipara)
     end do
  end do

  ! Reading input for bond angle parameters about y2
  rewind(lunpara)
  
  call ukoto_uiread2(lunpara, lunout, 'bond_angle_y2   ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "bond_angle_y2" in the flexible_local.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     read (ctmp00, *) amino1, &
          param(1), param(2), param(3), param(4), param(5), &
          &param(6), param(7), param(8), param(9), param(10)
     aaid(1) = ifunc_seq2id(amino1)
     
     do ipara = 1, npara_ang
        inflp%ang_para_y2(aaid(1), ipara) = param(ipara)
     end do
  end do

#ifdef MPI_PAR
  end if
  call MPI_Bcast (inflp, inflp%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_mapara_flp

