! setp_mapara_ele
!> @brief Read parameters from electrostatic.para file. The parameters are &
!>        used for electrostatic interaction

subroutine setp_mapara_ele()
  
  use const_maxsize
  use const_index
  use const_physical
  use var_inp, only : infile, outfile
  use var_setp, only : inele
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ----------------------------------------------------------------------
  ! intent(out) :: inpara

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: i, aaid
  integer :: lunout, lunpara
  integer :: iline, nlines, iequa, nequat
  real(PREC) :: x

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
  lunpara = infile%para_ele

  ! -------------------------------------------------------------------
  inele%i_diele  = 0   ! default
  inele%i_charge = 0   ! default
  inele%i_function_form = 0   ! default
  inele%cutoff_ele     = INVALID_VALUE
  inele%ionic_strength = INVALID_VALUE
  inele%diele_water    = INVALID_VALUE
  inele%coef_charge_type(:) = INVALID_VALUE
  inele%length_per_unit(:) = 0.0  ! default
  inele%dna2_phos_pro_charge = -1.0   ! default


  ! -------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'para_cafemol_ele', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "para_cafemol_ele" in the electrostatic%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat
        cvalue = 'cutoff_ele'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inele%cutoff_ele, cvalue)

        cvalue = 'ionic_strength'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inele%ionic_strength, cvalue)
        
        cvalue = 'diele_water'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inele%diele_water, cvalue)
        
        cvalue = 'i_diele'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inele%i_diele, cvalue)

        cvalue = 'dna2_phos_pro_charge'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inele%dna2_phos_pro_charge, cvalue)

        cvalue = 'i_charge'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inele%i_charge, cvalue)

        cvalue = 'i_function_form'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inele%i_function_form, cvalue)
     end do ! iequa

     if(ctmp00(1:11) == 'CHARGE_TYPE') then
        char3 = '   '
        read (ctmp00, *) char50, char3, x
        write (lunout, '(5a, g10.3)') '---reading CHARGE_TYPE: ', trim(char50), ' ', char3, ' ', x

        aaid = ifunc_seq2id(char3)
        if(aaid < 1 .or. aaid > SEQID%CL) then
           write(error_message,*) 'Error: invalid charge type is specified aaid=', aaid
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        inele%coef_charge_type(aaid) = x

     end if

     if(ctmp00(1:15) == 'LENGTH_PER_UNIT') then
        char3 = '   '
        read (ctmp00, *) char50, char3, x
        write (lunout, '(5a, g10.3)') '---reading LENGTH_PER_UNIT: ', trim(char50), ' ', char3, ' ', x

        if (char3 == 'OTH') then
            aaid = 21
        else if (char3 == 'P  ') then
            aaid = 22
        else
            aaid = ifunc_seq2id(char3)
        end if
        inele%length_per_unit(aaid) = x
     end if
  end do ! iline

  ! -------------------------------------------------------------------
  ! check
  if (inele%ionic_strength > INVALID_JUDGE) then
     error_message = 'Error: invalid value for inele%ionic_strength'
     call util_error(ERROR%STOP_ALL, error_message)
  endif

  if (inele%cutoff_ele > INVALID_JUDGE) then
     error_message = 'Error: invalid value for inele%cutoff_ele'
     call util_error(ERROR%STOP_ALL, error_message)
  endif

  if (inele%ionic_strength < ZERO_JUDGE) then
     error_message = 'Error: ionic_strength is too small'
     call util_error(ERROR%STOP_ALL, error_message)
  endif

  if(inele%i_diele == 0) then
!     write (lunout, *) 'using constant value for dielectric constant'
     if (inele%diele_water > INVALID_JUDGE) then
        error_message = 'Error: invalid value for inele%diele_water'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  else if (inele%i_diele == 1) then
!     write (lunout, *) 'using dielectric constant as function of temperature and concentration of cation'
  else if (inele%i_diele == 2) then
     continue
  else
     error_message = 'Error: invalid value for inele%i_diele'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  if(inele%i_charge == 0) then
     continue  ! default
  else if (inele%i_charge == 1) then
     !write (lunout, *) "using charge values based on Manning's theory"
  else
     error_message = 'Error: invalid value for inele%i_charge'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  if(inele%i_function_form == 0) then
     continue  ! default
  else if (inele%i_function_form == 1) then
     !write (lunout, *) "Coulomb potential is enabled."
  else
     error_message = 'Error: invalid value for inele%i_function_form'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do i = 1, SEQID%CL
     if(inele%coef_charge_type(i) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for coef_charge_type'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
  end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast (inele,inele%sz,MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_mapara_ele
