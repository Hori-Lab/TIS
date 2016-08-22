! setp_mapara_aicg
!> @brief Subroutine for reading parameters for AICG

! ***********<subrotine for the common about parameter>*************
! This subrouine is for reading the paramter sets.
subroutine setp_mapara_aicg()
  
  use const_maxsize
  use const_index
  use var_io,  only : infile, outfile
  use var_setp, only : inaicg, inaicg2
  use mpiconst

  implicit none

  integer :: lunout, lunpara
  integer :: iline, nlines, iequa, nequat
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  lunout = outfile%data
  lunpara = infile%para_aicg_gen

  ! -------------------------------------------------------------------
  ! aicg
  inaicg%cbd_aicg = -1.0
  inaicg%cba_aicg_G = -1.0
  inaicg%cba_aicg_H = -1.0
  inaicg%cba_aicg_E = -1.0
  inaicg%cba_aicg_T = -1.0
  inaicg%cba_aicg_C = -1.0
  inaicg%cdih_aicg_G = -1.0
  inaicg%cdih_aicg_H = -1.0
  inaicg%cdih_aicg_E = -1.0
  inaicg%cdih_aicg_T = -1.0
  inaicg%cdih_aicg_C = -1.0
  inaicg%ave_caicg = -1.0
  inaicg%gen_caicg = -1.0
  inaicg%ecut_low = -1.0
  inaicg%ecut_up = -1.0
  inaicg%iflag_scale = -1
  ! aicg2
  inaicg2%cbd_aicg2 = -1.0
  inaicg2%ave_caicg2_loc = -1.0
  inaicg2%gen_caicg2_loc = -1.0
  inaicg2%wid_aicg13 = -1.0
  inaicg2%wid_aicg14 = -1.0
  inaicg2%wid_dih = -1.0
  inaicg2%ave_caicg2plus_13 = -1.0
  inaicg2%gen_caicg2plus_13 = -1.0
  inaicg2%ave_caicg2plus_14 = -1.0
  inaicg2%gen_caicg2plus_14 = -1.0
  inaicg2%ave_caicg2plus_nloc = -1.0
  inaicg2%gen_caicg2plus_nloc = -1.0
  inaicg2%ave_caicg2_nloc = -1.0
  inaicg2%gen_caicg2_nloc = -1.0
  inaicg2%ecut_low_aicg2 = -1.0
  inaicg2%ecut_up_aicg2 = -1.0
  inaicg2%iflag_scale_aicg2 = -1

  ! ------------------------------------------------------------------- 

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'aicg_para        ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "aicg_para       " field in the infile%para_aicg_gen file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat
        
        cvalue = 'cbd_aicg'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%cbd_aicg, cvalue)
        
        cvalue = 'cba_aicg_G'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%cba_aicg_G, cvalue)
        
        cvalue = 'cba_aicg_H'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%cba_aicg_H, cvalue)
        
        cvalue = 'cba_aicg_E'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%cba_aicg_E, cvalue)
        
        cvalue = 'cba_aicg_T'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%cba_aicg_T, cvalue)
        
        cvalue = 'cba_aicg_C'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%cba_aicg_C, cvalue)
        
        cvalue = 'cdih_aicg_G'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%cdih_aicg_G, cvalue)
        
        cvalue = 'cdih_aicg_H'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%cdih_aicg_H, cvalue)
        
        cvalue = 'cdih_aicg_E'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%cdih_aicg_E, cvalue)
        
        cvalue = 'cdih_aicg_T'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%cdih_aicg_T, cvalue)
        
        cvalue = 'cdih_aicg_C'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%cdih_aicg_C, cvalue)
        
        cvalue = 'ave_caicg'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%ave_caicg, cvalue)
        
        cvalue = 'gen_caicg'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%gen_caicg, cvalue)
        
        cvalue = 'ecut_low'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%ecut_low, cvalue)
        
        cvalue = 'ecut_up'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg%ecut_up, cvalue)
        
        cvalue = 'iflag_scale'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inaicg%iflag_scale, cvalue)
     end do
  end do

  ! -------------------------------------------------------------------
  if(inaicg%cbd_aicg < 0.0) then
     error_message = 'Error: invalid value for cbd_aicg'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%cba_aicg_G < 0.0) then
     error_message = 'Error: invalid value for cba_aicg_G'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%cba_aicg_H < 0.0) then
     error_message = 'Error: invalid value for cba_aicg_H'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%cba_aicg_E < 0.0) then
     error_message = 'Error: invalid value for cba_aicg_E'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%cba_aicg_T < 0.0) then
     error_message = 'Error: invalid value for cba_aicg_T'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%cba_aicg_C < 0.0) then
     error_message = 'Error: invalid value for cba_aicg_C'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%cdih_aicg_G < 0.0) then
     error_message = 'Error: invalid value for cdih_aicg_G'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%cdih_aicg_H < 0.0) then
     error_message = 'Error: invalid value for cdih_aicg_H'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%cdih_aicg_E < 0.0) then
     error_message = 'Error: invalid value for cdih_aicg_E'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%cdih_aicg_T < 0.0) then
     error_message = 'Error: invalid value for cdih_aicg_T'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%cba_aicg_C < 0.0) then
     error_message = 'Error: invalid value for cba_aicg_C'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%ave_caicg < 0.0) then
     error_message = 'Error: invalid value for ave_caicg'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%gen_caicg < 0.0) then
     error_message = 'Error: invalid value for gen_caicg'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%ecut_low < 0.0) then
     error_message = 'Error: invalid value for ecut_low'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%ecut_up < 0.0) then
     error_message = 'Error: invalid value for ecut_up'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg%iflag_scale < 0 .or. inaicg%iflag_scale > 1) then
     error_message = 'Error: invalid value for iflag_scale'
     call util_error(ERROR%STOP_ALL, error_message)

  end if

  ! aicg2
  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'aicg2_para        ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "aicg2_para       " field in the aicg_generic para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat
        
        cvalue = 'cbd_aicg2'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%cbd_aicg2, cvalue)
        
        cvalue = 'wid_aicg13'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%wid_aicg13, cvalue)
        
        cvalue = 'wid_aicg14'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%wid_aicg14, cvalue)

        cvalue = 'ave_caicg2_loc'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%ave_caicg2_loc, cvalue)
        
        cvalue = 'gen_caicg2_loc'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%gen_caicg2_loc, cvalue)
        
        cvalue = 'ave_caicg2_nloc'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%ave_caicg2_nloc, cvalue)

        cvalue = 'gen_caicg2_nloc'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%gen_caicg2_nloc, cvalue)

        cvalue = 'wid_dih'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%wid_dih, cvalue)

        cvalue = 'ave_caicg2plus_13'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%ave_caicg2plus_13, cvalue)
        
        cvalue = 'gen_caicg2plus_13'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%gen_caicg2plus_13, cvalue)
        
        cvalue = 'ave_caicg2plus_14'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%ave_caicg2plus_14, cvalue)
  
        cvalue = 'gen_caicg2plus_14'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%gen_caicg2plus_14, cvalue)

        cvalue = 'ave_caicg2plus_nloc'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%ave_caicg2plus_nloc, cvalue)

        cvalue = 'gen_caicg2plus_nloc'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%gen_caicg2plus_nloc, cvalue)
        
        cvalue = 'ecut_low_aicg2'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%ecut_low_aicg2, cvalue)

        cvalue = 'ecut_up_aicg2'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inaicg2%ecut_up_aicg2, cvalue)
        
        cvalue = 'iflag_scale_aicg2'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inaicg2%iflag_scale_aicg2, cvalue)
     end do
  end do

  ! -------------------------------------------------------------------
  if(inaicg2%cbd_aicg2 < 0.0) then
     error_message = 'Error: invalid value for cbd_aicg2'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%wid_aicg13 < 0.0) then
     error_message = 'Error: invalid value for wid_aicg13'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%wid_aicg14 < 0.0) then
     error_message = 'Error: invalid value for wid_aicg14'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%wid_dih < 0.0) then
     error_message = 'Error: invalid value for wid_dih'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%ave_caicg2_loc < 0.0) then
     error_message = 'Error: invalid value for ave_caicg2_loc'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%ave_caicg2plus_13 < 0.0) then
     error_message = 'Error: invalid value for ave_caicg2plus_13'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%ave_caicg2plus_14 < 0.0) then
     error_message = 'Error: invalid value for ave_caicg2plus_14'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%ave_caicg2_nloc < 0.0) then
     error_message = 'Error: invalid value for ave_caicg2_nloc'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%ave_caicg2plus_nloc < 0.0) then
     error_message = 'Error: invalid value for ave_caicg2plus_nloc'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%gen_caicg2_loc < 0.0) then
     error_message = 'Error: invalid value for gen_caicg2_loc'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%gen_caicg2plus_13 < 0.0) then
     error_message = 'Error: invalid value for gen_caicg2plus_13'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%gen_caicg2plus_14 < 0.0) then
     error_message = 'Error: invalid value for gen_caicg2plus_14'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%gen_caicg2_nloc < 0.0) then
     error_message = 'Error: invalid value for gen_caicg2_nloc'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%gen_caicg2plus_nloc < 0.0) then
     error_message = 'Error: invalid value for gen_caicg2plus_nloc'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%ecut_low_aicg2 < 0.0) then
     error_message = 'Error: invalid value for ecut_low_aicg2'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%ecut_up_aicg2 < 0.0) then
     error_message = 'Error: invalid value for ecut_up_aicg2'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inaicg2%iflag_scale_aicg2 < 0 .or. inaicg2%iflag_scale_aicg2 > 1) then
     error_message = 'Error: invalid value for iflag_scale_aicg2'
     call util_error(ERROR%STOP_ALL, error_message)

  end if
#ifdef MPI_PAR
  end if

  call MPI_Bcast (inaicg, inaicg%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (inaicg2,inaicg2%sz,MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_mapara_aicg
