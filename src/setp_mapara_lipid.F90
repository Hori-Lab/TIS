! setp_mapara_lipid
!> @brief This subroutine is to read the parameters related to lipid simulation from "para_lip" file. &
!>        All the information are stored into "inlip" struct.

! ***********<subrotine for the common about parameter>*************
! This subrouine is for reading the paramter sets.
subroutine setp_mapara_lipid()
  
  use const_maxsize
  use const_index
  use var_inp,  only : outfile, infile
  use var_setp, only : inlip
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
  lunpara = infile%para_lip

  ! -------------------------------------------------------------------
  inlip%energy_unit_lipid = -1.0
  inlip%num_lip_total     = -1
  inlip%num_lip_core      = -1
  inlip%num_lip_int       = -1
  inlip%num_lip_tail      = -1
  inlip%sigma_lipid       = -1.0
  inlip%cbd_lipid         = -1.0
  inlip%cba_lipid         = -1.0
  inlip%ccore         = -1.0
  inlip%cutoff_core       = -1.0
  inlip%ctail         = -1.0
  inlip%cutoff_tail       = -1.0
  inlip%cint          = -1.0
  inlip%cutoff_int        = -1.0

  ! -------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'para_cafemol_lip', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "para_cafemol_lip" field in the lipid%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat
        cvalue = 'energy_unit_lipid'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%energy_unit_lipid, cvalue)
        
        cvalue = 'num_lip_total'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inlip%num_lip_total, cvalue)
        
        cvalue = 'num_lip_core'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inlip%num_lip_core, cvalue)
        
        cvalue = 'num_lip_int'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inlip%num_lip_int, cvalue)
        
        cvalue = 'num_lip_tail'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inlip%num_lip_tail, cvalue)
        
        cvalue = 'sigma_lipid'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%sigma_lipid, cvalue)
        
        cvalue = 'cbd_lipid'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cbd_lipid, cvalue)
        
        cvalue = 'cba_lipid'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cba_lipid, cvalue)
        
        cvalue = 'ccore'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%ccore, cvalue)
        
        cvalue = 'cutoff_core'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cutoff_core, cvalue)
        
        cvalue = 'ctail'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%ctail, cvalue)
        
        cvalue = 'cutoff_tail'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cutoff_tail, cvalue)
        
        cvalue = 'cint'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cint, cvalue)
        
        cvalue = 'cutoff_int'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cutoff_int, cvalue)
     end do
  end do
  

  ! -------------------------------------------------------------------
  if(inlip%energy_unit_lipid < 0.0) then
     error_message = 'Error: invalid value for energy_unit_lipid'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%num_lip_total < 0) then
     error_message = 'Error: invalid value for num_lip_total'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%num_lip_core < 0) then
     error_message = 'Error: invalid value for num_lip_core'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%num_lip_int < 0) then
     error_message = 'Error: invalid value for num_lip_int'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%num_lip_tail < 0) then
     error_message = 'Error: invalid value for num_lip_tail'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%sigma_lipid < 0.0) then
     error_message = 'Error: invalid value for sigma_lipid'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%cbd_lipid < 0.0) then
     error_message = 'Error: invalid value for cbd_lipid'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%cba_lipid < 0.0) then
     error_message = 'Error: invalid value for cba_lipid'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%ccore < 0.0) then
     error_message = 'Error: invalid value for ccore'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%cutoff_core < 0.0) then
     error_message = 'Error: invalid value for cutoff_core'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%ctail < 0.0) then
     error_message = 'Error: invalid value for ctail'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%cutoff_tail < 0.0) then
     error_message = 'Error: invalid value for cutoff_tail'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%cint < 0.0) then
     error_message = 'Error: invalid value for cint'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inlip%cutoff_int < 0.0) then
     error_message = 'Error: invalid value for cutoff_int'
     call util_error(ERROR%STOP_ALL, error_message)
     
  end if

  ! -------------------------------------------------------------------
  if(inlip%num_lip_total /= inlip%num_lip_core + inlip%num_lip_int + inlip%num_lip_tail) then
     error_message = 'Error: num_lip_total should be equal num_lip_core + num_lip_int + num_lip_tail'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast (inlip, inlip%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_mapara_lipid
