! setp_mapara_hp
!> @brief Read parameters from hydrophobic.para file. The parameters are &
!>        used for hydrophobic interaction

subroutine setp_mapara_hp()
  
  use const_maxsize
  use const_index
  use const_physical
  use var_io, only : infile, outfile
  use var_setp, only : inhp
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ----------------------------------------------------------------------
  ! intent(out) :: inpara

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: i, aaid, j
  integer :: lunout, lunpara
  integer :: iline, nlines, iequa, nequat
  integer :: num(21)
  integer :: ncoor
  real(PREC) :: ncoormax
  real(PREC) :: chp

  character(3) :: char3
  character(50) :: char50
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  integer :: ifunc_seq2id

  ! -------------------------------------------------------------------
  lunout = outfile%data
  lunpara = infile%para_hp

  ! -------------------------------------------------------------------
  inhp%rho_min_hp = -1.0
  inhp%coef_rho_hp = -1.0
  inhp%coef_hp = -1.0
  inhp%ncoor_para_hp(1:21) = -1
  inhp%ncoormax_para_hp(1:21) = -1
  inhp%coefaa_para_hp(1:21) = -1.0
  inhp%cutoffdmin_para_hp(1:21, 1:21) = -1.0
  inhp%cutoffdmax_para_hp(1:21, 1:21) = -1.0


  ! -------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'para_cafemol_hp ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "para_cafemol_hp" in the hydrophobic%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat
        cvalue = 'coef_hp'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inhp%coef_hp, cvalue)

        cvalue = 'coef_rho_hp'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inhp%coef_rho_hp, cvalue)
        
        cvalue = 'rho_min_hp'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inhp%rho_min_hp, cvalue)
     end do ! iequa

     if(ctmp00(1:8) == 'HPE_TYPE') then
        read (ctmp00, *) char50, char3, chp, ncoor, ncoormax
        write (lunout, '(3a, f6.2,i5,f6.2)') '---reading HPE_TYPE: ', trim(char50), char3, chp, ncoor, ncoormax
        if (char3 == 'OTH') then
            aaid = 21
        else
            aaid = ifunc_seq2id(char3)
        end if
        inhp%coefaa_para_hp(aaid) = chp
        inhp%ncoor_para_hp(aaid) = ncoor 
        inhp%ncoormax_para_hp(aaid) = ncoormax
     end if
  end do ! iline

  ! ---------------------------------------------------------------------------------
  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'cutoff_dmin_hp  ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "cutoff_dmin_hp" in the hydrophobic%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  num(1:21) = 0
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     if(ctmp00(1:10) == 'HPE_CUTOFF') then
        read (ctmp00, *) char50, char3 
        if (char3 == 'OTH') then
            aaid = 21
        else
            aaid = ifunc_seq2id(char3)
        end if
        read (ctmp00, *) char50, char3, &
          (inhp%cutoffdmin_para_hp(aaid, j+num(aaid)), j=1, 7)
        write (lunout, '(3a, 7g10.3)') '---reading HPE_CUTOFF: ', &
             trim(char50), char3, &
          (inhp%cutoffdmin_para_hp(aaid, j+num(aaid)), j=1, 7)
        num(aaid) = num(aaid) + 7
     end if
  end do

  ! ---------------------------------------------------------------------------------
  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'cutoff_dmax_hp  ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "cutoff_dmax_hp" in the hydrophobic%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  num(1:21) = 0
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     if(ctmp00(1:10) == 'HPE_CUTOFF') then
        read (ctmp00, *) char50, char3 
        if (char3 == 'OTH') then
            aaid = 21
        else
            aaid = ifunc_seq2id(char3)
        end if
        read (ctmp00, *) char50, char3, &
             (inhp%cutoffdmax_para_hp(aaid, j+num(aaid)), j=1, 7 ) 
        write (lunout, '(3a, 7g10.3)') '---reading HPE_CUTOFF: ', &
             trim(char50), char3, &
             (inhp%cutoffdmax_para_hp(aaid, j+num(aaid)), j=1, 7)
        num(aaid) = num(aaid) + 7 
     end if
  end do

  ! -------------------------------------------------------------------
  if(inhp%coef_hp < 0.0) then
     error_message = 'Error: invalid value for coef_hp'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inhp%coef_rho_hp < 0.0) then
     error_message = 'Error: invalid value for coef_rho_hp'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inhp%rho_min_hp < 0.0) then
     error_message = 'Error: invalid value for rho_min_hp'
     call util_error(ERROR%STOP_ALL, error_message)
     
  end if
  do i = 1, 21 
     if(inhp%ncoor_para_hp(i) < 0) then
        error_message = 'Error: invalid value for ncoor_para_hp'
        call util_error(ERROR%STOP_ALL, error_message)
     else if(inhp%ncoormax_para_hp(i) < 0) then
        error_message = 'Error: invalid value for ncoormax_para_hp'
        call util_error(ERROR%STOP_ALL, error_message)
!     else if(inhp%coefaa_para_hp(i) < 0.0) then
!        error_message = 'Error: invalid value for coefaa_para_hp'
!        call util_error(ERROR%STOP_ALL, error_message)
     else
        do j = 1, 21 
           if(inhp%cutoffdmin_para_hp(i,j) < 0.0) then
              error_message = 'Error: invalid value for cutoffdmin_para_hp'
              call util_error(ERROR%STOP_ALL, error_message)
           else if(inhp%cutoffdmax_para_hp(i,j) < 0.0) then
              error_message = 'Error: invalid value for cutoffdmax_para_hp'
              call util_error(ERROR%STOP_ALL, error_message)
           end if
        end do
     end if
  end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast (inhp,inhp%sz,MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_mapara_hp
