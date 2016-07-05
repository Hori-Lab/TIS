!setp_mapara_pro
!> @brief Reads parameters related to protein simulation from a parameter &
!>        file for protein.                                               &
!>        All the information are store into "inpro" struct.

subroutine setp_mapara_pro(lunpara, lunout)
  
  use const_maxsize
  use const_index
  use var_setp, only : inpro
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ----------------------------------------------------------------------
  ! intent(out) :: inpro
  integer, intent(in) :: lunpara
  integer, intent(in) :: lunout

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: iline, nlines, iequa, nequat
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  inpro%energy_unit_protein = -1.0
  inpro%cbd                 = -1.0
  inpro%cba                 = -1.0
  inpro%cdih_1              = -1.0
  inpro%cdih_3              = -1.0
  inpro%n_sep_nlocal        = -1
  inpro%n_sep_contact       = -1
  inpro%cutoff_go           = -1.0
  inpro%cutoff_exvol        = -1.0
  inpro%dfcontact           = -1.0
  inpro%cgo1210             = -1.0
  inpro%cdist_rep12         = -1.0 
  inpro%cdist_rep6          = -1.0 
  inpro%crep12              = -1.0
  inpro%crep6               = -1.0

  ! ------------------------------------------------------------------- 

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'para_cafemol_pro', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "para_cafemol_pro" field in the protein%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat
        
        cvalue = 'energy_unit_protein'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%energy_unit_protein, cvalue)
        
        cvalue = 'cbd'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cbd, cvalue)
        
        cvalue = 'cba'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cba, cvalue)
        
        cvalue = 'cdih_1'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cdih_1, cvalue)
        
        cvalue = 'cdih_3'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cdih_3, cvalue)
        
        cvalue = 'n_sep_nlocal'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inpro%n_sep_nlocal, cvalue)
        
        cvalue = 'n_sep_contact'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inpro%n_sep_contact, cvalue)
        
        cvalue = 'cutoff_go'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cutoff_go, cvalue)

        cvalue = 'cutoff_exvol'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cutoff_exvol, cvalue)

        cvalue = 'dfcontact'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%dfcontact, cvalue)
        
        cvalue = 'cgo1210'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cgo1210, cvalue)
        
        cvalue = 'cdist_rep12'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cdist_rep12, cvalue)
        
        cvalue = 'crep12'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%crep12, cvalue)

        cvalue = 'cdist_rep6'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cdist_rep6, cvalue)
        
        cvalue = 'crep6'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%crep6, cvalue)
     end do
  end do

  ! -------------------------------------------------------------------
  if(inpro%energy_unit_protein < 0.0) then
     error_message = 'Error: invalid value for energy_unit_protein'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%cbd < 0.0) then
     error_message = 'Error: invalid value for cbd'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%cba < 0.0) then
     error_message = 'Error: invalid value for cba'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%cdih_1 < 0.0) then
     error_message = 'Error: invalid value for cdih_1'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%cdih_3 < 0.0) then
     error_message = 'Error: invalid value for cdih_3'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%n_sep_nlocal < 0) then
     error_message = 'Error: invalid value for n_sep_nlocal'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%n_sep_contact < 0) then
     error_message = 'Error: invalid value for n_sep_contat'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%cutoff_go < 0.0) then
     error_message = 'Error: invalid value for cutoff_go'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%cutoff_exvol < 0.0) then
     error_message = 'Error: invalid value for cutoff_exvol'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%dfcontact < 0.0) then
     error_message = 'Error: invalid value for dfcontact'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%cgo1210 < 0.0) then
     error_message = 'Error: invalid value for cgo1210'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%cdist_rep12 < 0.0) then
     error_message = 'Error: invalid value for cdist_rep12'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%crep12 < 0.0) then
     error_message = 'Error: invalid value for crep12'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%cdist_rep6 < 0.0) then
     error_message = 'Error: invalid value for cdist_rep6'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpro%crep6 < 0.0) then
     error_message = 'Error: invalid value for crep6'
     call util_error(ERROR%STOP_ALL, error_message)

  end if


  ! -------------------------------------------------------------------
  if(inpro%n_sep_contact < inpro%n_sep_nlocal) then
     error_message = 'Error: n_sep_contact must be larger than n_sep_nlocal'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast (inpro, inpro%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_mapara_pro
