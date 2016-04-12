!setp_mapara_dna
!> @brief Reads parameters related to DNA simulation from "para_dna" file. &
!>        All the information are store into "indna" struct.

subroutine setp_mapara_dna()
  
  use const_maxsize
  use const_index
  use var_inp, only : infile, outfile
  use var_setp, only : indna
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
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  lunout = outfile%data
  lunpara = infile%para_dna

  ! -------------------------------------------------------------------
  indna%energy_unit_dna   = -1.0
  indna%cbd_dna           = -1.0
  indna%cbd2_dna          = -1.0
  indna%cba_dna           = -1.0
  indna%cdih_1_dna        = -1.0
  indna%cdih_3_dna        = -1.0
  indna%dfcontact_dna     = -1.0
  indna%renorm_dist_stack = -1.0
  indna%cstack        = -1.0
  indna%cutoff_stack      = -1.0
  indna%cbp_at        = -1.0
  indna%cdist_bp_at       = -1.0
  indna%cbp_gc        = -1.0
  indna%cdist_bp_gc       = -1.0
  indna%cutoff_bp         = -1.0
  indna%cmbp          = -1.0
  indna%cdist_mbp         = -1.0
  indna%cutoff_mbp        = -1.0
  indna%cexv_dna      = -1.0
  indna%cdist_exv_dna     = -1.0
  indna%cutoff_exv_dna    = -1.0
  indna%csolvmax_dna  = -1.0
  indna%cdist_solv_dna    = -1.0
  indna%cralpha_solv_dna  = -1.0
  indna%cutoff_solv_dna   = -1.0


  ! -------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'para_cafemol_dna', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "para_cafemol_dna" in the dna%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat
        cvalue = 'energy_unit_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%energy_unit_dna, cvalue)
        
        cvalue = 'cbd_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cbd_dna, cvalue)
        
        cvalue = 'cbd2_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cbd2_dna, cvalue)
        
        cvalue = 'cba_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cba_dna, cvalue)
        
        cvalue = 'cdih_1_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdih_1_dna, cvalue)
        
        cvalue = 'cdih_3_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdih_3_dna, cvalue)
        
        cvalue = 'dfcontact_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%dfcontact_dna, cvalue)
        
        cvalue = 'renorm_dist_stack'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%renorm_dist_stack, cvalue)
        
        cvalue = 'cstack'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cstack, cvalue)
        
        cvalue = 'cutoff_stack'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cutoff_stack, cvalue)
        
        cvalue = 'cbp_at'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cbp_at, cvalue)
        
        cvalue = 'cdist_bp_at'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdist_bp_at, cvalue)
        
        cvalue = 'cbp_gc'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cbp_gc, cvalue)
        
        cvalue = 'cdist_bp_gc'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdist_bp_gc, cvalue)
        
        cvalue = 'cutoff_bp'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cutoff_bp, cvalue)
        
        cvalue = 'cmbp'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cmbp, cvalue)
        
        cvalue = 'cdist_mbp'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdist_mbp, cvalue)
        
        cvalue = 'cutoff_mbp'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cutoff_mbp, cvalue)
        
        cvalue = 'cexv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cexv_dna, cvalue)
        
        cvalue = 'cdist_exv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdist_exv_dna, cvalue)
        
        cvalue = 'cutoff_exv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cutoff_exv_dna, cvalue)
        
        cvalue = 'csolvmax_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%csolvmax_dna, cvalue)
        
        cvalue = 'cdist_solv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdist_solv_dna, cvalue)
        
        cvalue = 'cralpha_solv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cralpha_solv_dna, cvalue)
        
        cvalue = 'cutoff_solv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cutoff_solv_dna, cvalue)
     end do
  end do

  ! -------------------------------------------------------------------
  if(indna%energy_unit_dna < 0.0) then
     error_message = 'Error: invalid value for energy_unit_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cbd_dna < 0.0) then
     error_message = 'Error: invalid value for cbd_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cbd2_dna < 0.0) then
     error_message = 'Error: invalid value for cbd2_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cba_dna < 0.0) then
     error_message = 'Error: invalid value for cba_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cdih_1_dna < 0.0) then
     error_message = 'Error: invalid value for cdih_1_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cdih_3_dna < 0.0) then
     error_message = 'Error: invalid value for cdih_3_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%dfcontact_dna < 0.0) then
     error_message = 'Error: invalid value for dfcontact_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%renorm_dist_stack < 0.0) then
     error_message = 'Error: invalid value for renorm_dist_stack'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cstack < 0.0) then
     error_message = 'Error: invalid value for cstack'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cutoff_stack < 0.0) then
     error_message = 'Error: invalid value for cutoff_stack'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cbp_at < 0.0) then
     error_message = 'Error: invalid value for cbp_at'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cdist_bp_at < 0.0) then
     error_message = 'Error: invalid value for cdist_bp_at'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cbp_gc < 0.0) then
     error_message = 'Error: invalid value for cbp_gc'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cdist_bp_gc < 0.0) then
     error_message = 'Error: invalid value for cdist_bp_gc'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cutoff_bp < 0.0) then
     error_message = 'Error: invalid value for cutoff_bp'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cmbp < 0.0) then
     error_message = 'Error: invalid value for cmbp'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cdist_mbp < 0.0) then
     error_message = 'Error: invalid value for cdist_mbp'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cutoff_mbp < 0.0) then
     error_message = 'Error: invalid value for cutoff_mbp'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cexv_dna < 0.0) then
     error_message = 'Error: invalid value for cexv_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cdist_exv_dna < 0.0) then
     error_message = 'Error: invalid value for cdist_exv_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cutoff_exv_dna < 0.0) then
     error_message = 'Error: invalid value for cutoff_exv_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%csolvmax_dna < 0.0) then
     error_message = 'Error: invalid value for csolvmax_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cdist_solv_dna < 0.0) then
     error_message = 'Error: invalid value for cdist_solv_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cralpha_solv_dna < 0.0) then
     error_message = 'Error: invalid value for cralpha_solv_dna'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(indna%cutoff_solv_dna < 0.0) then
     error_message = 'Error: invalid value for cutoff_solv_dna'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
#ifdef MPI_PAR
  end if

  call MPI_Bcast (indna, indna%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_mapara_dna
