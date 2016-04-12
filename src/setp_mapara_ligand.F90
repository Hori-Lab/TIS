! setp_mapara_ligand
!> @brief Read parameters from ligand.para file. The parameters are &
!>        used for explicit ligand model

subroutine setp_mapara_ligand()
  
  use const_maxsize
  use const_index
  use const_physical
  use var_inp, only : infile, outfile
  use var_setp, only : inligand
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ----------------------------------------------------------------------
  ! intent(out) :: inpara

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: lunout, lunpara
  integer :: iline, nlines, iequa, nequat
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  lunout = outfile%data
  lunpara = infile%para_lig

  ! -------------------------------------------------------------------'
  inligand%energy_unit = -1.0
  inligand%cbd = -1.0
  inligand%cba = -1.0
  inligand%cdih = -1.0
  inligand%crep12 = -1.0
  inligand%cdist_rep12_lpro = -1.0
  inligand%cdist_rep12_llig = -1.0
  inligand%cutoff_exvol = -1.0

  ! -------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'para_cafemol_lig', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "para_cafemol_lig" in the ligand%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat
        cvalue = 'energy_unit_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%energy_unit, cvalue)
        
        cvalue = 'cbd_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cbd, cvalue)
        
        cvalue = 'cba_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cba, cvalue)
        
        cvalue = 'cdih_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cdih, cvalue)
        
        cvalue = 'crep12_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%crep12, cvalue)
        
        cvalue = 'cdist_rep12_lpro'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cdist_rep12_lpro, cvalue)
        
        cvalue = 'cdist_rep12_llig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cdist_rep12_llig, cvalue)
        
        cvalue = 'cutoff_exvol_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cutoff_exvol, cvalue)
     end do ! iequa

  end do ! iline

  ! -------------------------------------------------------------------
  if(inligand%energy_unit < 0.0) then
     error_message = 'Error: invalid value for energy_unit_lig'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inligand%cbd < 0.0) then
     error_message = 'Error: invalid value for cbd_lig'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inligand%cba < 0.0) then
     error_message = 'Error: invalid value for cba_lig'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inligand%cdih < 0.0) then
     error_message = 'Error: invalid value for cdih_lig'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inligand%crep12 < 0.0) then
     error_message = 'Error: invalid value for crep12_lig'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inligand%cdist_rep12_lpro < 0.0) then
     error_message = 'Error: invalid value for cdist_rep12_lpro'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inligand%cdist_rep12_llig < 0.0) then
     error_message = 'Error: invalid value for cdist_rep12_llig'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inligand%cutoff_exvol < 0.0) then
     error_message = 'Error: invalid value for cutoff_exvol_lig'
     call util_error(ERROR%STOP_ALL, error_message)
  endif
   
#ifdef MPI_PAR
  end if

  call MPI_Bcast (inligand,inligand%sz,MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_mapara_ligand
