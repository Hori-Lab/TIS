! setp_para_lipid
!> @brief Reads "initial_lipid" block in input-file and stores the data into  &
!>        "inlip" struct.

! ******************************************************************
subroutine setp_para_lipid()
  
  use const_maxsize
  use const_index
  use var_inp, only : infile, outfile
  use var_setp, only: inlip

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  !------------------------------------------------------------------
  ! local variables
  integer :: ilayer
  integer :: icol
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
  inlip%nmp_transverse_lipid = -1
  inlip%nmp_longitudinal_lipid = -1
  inlip%nlayer_lipid = -1
  inlip%grid_size_lipid = -1
  inlip%z_coord_lipid(1:MXLAYER) = -1000000.0

  ! ---------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'initial_lipid   ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "initial_lipid" field in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
         
     do iequa = 1, nequat
        ctmp02 = csides(1, iequa)

        cvalue = 'nmp_transverse_lipid'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inlip%nmp_transverse_lipid, cvalue)

        cvalue = 'nmp_longitudinal_lipid'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inlip%nmp_longitudinal_lipid, cvalue)

        cvalue = 'nlayer_lipid'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inlip%nlayer_lipid, cvalue)

        cvalue = 'grid_size_lipid'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%grid_size_lipid, cvalue)

        cvalue = 'z_coord_lipid'
        if(ctmp02(1:13) == cvalue) then
           do icol = 15, CARRAY_MXCOLM
              if(ctmp02(icol:icol) == ')') exit
           end do
           read (ctmp02(15:icol-1), *) ilayer
           read (csides(2, iequa), *) inlip%z_coord_lipid(ilayer)
        end if

     end do
  end do

  ! -----------------------------------------------------------------
  ! checking input variables
  if(inlip%nmp_transverse_lipid <= 0) then
     error_message = 'Error: invalid value for nmp_transverse_lipid'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inlip%nmp_longitudinal_lipid <= 0) then
     error_message = 'Error: invalid value for nmp_longitudinal_lipid'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inlip%nlayer_lipid <= 0) then
     error_message = 'Error: invalid value for nlayer_lipid'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inlip%nlayer_lipid > MXLAYER) then
     error_message = 'Error: should increase MXLAYER'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inlip%grid_size_lipid < 0.0) then
     error_message = 'Error: invalid value for grid_size_lipid'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do ilayer = 1, inlip%nlayer_lipid
     if(inlip%z_coord_lipid(ilayer) == -1000000.0) then
        error_message = 'Error: invalid value for z_coord_lipid'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
  end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast (inlip, inlip%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif


end subroutine setp_para_lipid
