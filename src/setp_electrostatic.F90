! setp_electrostatic
!> @brief This subroutine is to read and set the parameters for the electrostatic interaction.

subroutine setp_electrostatic()
  
  use const_maxsize
  use const_index
  use const_physical
  use var_io,  only : infile, outfile
  use var_struct, only : nunit_real, lunit2mp
  use var_setp, only : inele, inmisc, inpmf

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  integer :: ifunc_seq2id
  character(3) :: cfunc_id2seq

  ! -----------------------------------------------------------------------
  ! local variables
  integer :: i1, i2, k, aaid
  integer :: icol, isw, ifield, iulen
  integer :: iunit, junit, inunit(2), instate, iu, imp
  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  integer :: read_charge_change_imp
  integer :: i_charge_change
  logical :: icalc(MXUNIT)
  real(PREC) :: rcharge
  real(PREC) :: read_charge_change_value
  real(PREC) :: rvalue

  character(4) :: kfind
  character(12) :: char12
  character(13) :: header
  character(3) :: char3
  character(50) :: char50
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: ctmp02
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -----------------------------------------------------------------------
  iulen = 0
  luninp = infile%inp
  lunout = outfile%data

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

!  inele%i_diele  = 0   ! default
!  inele%cutoff_ele     = INVALID_VALUE
!  inele%ionic_strength = INVALID_VALUE
!  inele%diele_water    = INVALID_VALUE
  inele%i_calc_method = 0
  inele%n_charge_change = 0
  inele%i_function_form = -1
  inele%ewld_alpha = INVALID_VALUE
  inele%ewld_hmax = INVALID_VALUE
  inele%ewld_dipole = -1
  
  inele%flag_ele(1:MXMP) = .false.
  icalc(1:MXUNIT) = .false. 
  do iunit = 1, nunit_real
     do junit = iunit, nunit_real
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%ELE)) then
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
              inele%flag_ele(imp) = .true.
           end do
           do imp = lunit2mp(1, junit), lunit2mp(2, junit)
              inele%flag_ele(imp) = .true.
           end do
           icalc(iunit) = .true.
           icalc(junit) = .true.
        end if
     end do
  end do

  i_charge_change = 0

  ! -----------------------------------------------------------------------
  ! read from input file
  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'electrostatic   ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  if(kfind /= 'FIND') then
  !   error_message = 'Error: cannot find "electrostatic" field in the input file'
  !   call util_error(ERROR%STOP_ALL, error_message)

     ! Default parameters (para/electrostatic.para) are used.
     nlines = 0
  end if
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     ctmp00 = cwkinp(iline)

     do iequa = 1, nequat
        cvalue = 'cutoff_ele'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inele%cutoff_ele, cvalue)

        cvalue = 'ionic_strength'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inele%ionic_strength, cvalue)

        cvalue = 'conc_Mg'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inele%conc_Mg, cvalue)

        cvalue = 'diele_water'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inele%diele_water, cvalue)

        cvalue = 'i_diele'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inele%i_diele, cvalue)

        cvalue = 'i_charge'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inele%i_charge, cvalue)

        cvalue = 'i_function_form'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inele%i_function_form, cvalue)

        cvalue = 'i_calc_method'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inele%i_calc_method, cvalue)

        cvalue = 'ewld_alpha'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inele%ewld_alpha, cvalue)

        cvalue = 'ewld_hmax'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inele%ewld_hmax, cvalue)

        cvalue = 'ewld_dipole'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inele%ewld_dipole, cvalue)

        cvalue = 'i_semi'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inele%i_semi, cvalue)

        if(csides(1, iequa) == 'path_pmf_Mg_P') then
           inpmf%path(PMFTYPE%MG_P) = csides(2, iequa)
        end if

        cvalue = 'pmf_merge'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             rvalue, cvalue)
        inpmf%Rmerge(1:PMFTYPE%MAX) = rvalue
     end do

     if(ctmp00(1:11) == 'CHARGE_TYPE') then
        char3 = '   '
        read (ctmp00, *) char50, char3, rcharge
        write (lunout, '(5a, g10.3)') '---reading chage type: ', trim(char50), ' ', char3, ' ', rcharge
        aaid = ifunc_seq2id(char3)
        if(aaid < 1 .or. aaid > SEQID%CL) then
           write(error_message,*) 'Error: invalid charge type is specified aaid=', aaid
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        inele%coef_charge_type(aaid) = rcharge
        
     end if
  end do

  loop_line: do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     if(ctmp00(1:10) == 'CHARGE_DEL') then
        write (lunout, '(2a)') '---reading CHARGE_DEL: ', trim(ctmp00)
        isw = 0
        i1 = 0
        i2 = 0
        ifield = 0
        do k = 1, CARRAY_MXCOLM
           if(cwkinp(iline)(k:k) /= ' ' .and. &
                cwkinp(iline)(k:k) /= ',') then
              if(isw == 0) then
                 i1 = k
                 isw = 1
              end if
              i2 = k
           else if (isw == 1) then
              ifield = ifield + 1
              if(ifield == 2) then
                 read (ctmp00(i1:i2), *) char12
                 call util_unitstate(char12, inunit, instate)
                 iu = inunit(1)   !only consider real unit
                 if(iu <= 0 .or. iu > nunit_real) then
                    error_message = 'Error: invalid value for unit id in specifying charged sites'
                    call util_error(ERROR%STOP_ALL, error_message)
                 end if
                 iulen = lunit2mp(2, iu) - lunit2mp(1, iu) + 1
                 if(.not. icalc(iu)) then
                    write(error_message, '(a, i3, a)') 'Warning: unit ', iu,&
                    ' is not assigned electrostatic interaction '&
                    &//'but is specified the charged sites. Ignore'
                    call util_error(ERROR%WARN_ALL, error_message)
                    cycle loop_line 
                 end if
              else if(ifield == 3) then
                 read (ctmp00(i1:i2), *) char12
                 if(char12(1:1) /= 'u') then
                    error_message = 'Error: bad format for CHARGE_DEL line.  Forget the symbol u before particle id?'
                    call util_error(ERROR%STOP_ALL, error_message)
                 end if
              else if(ifield >= 4) then
                 read (ctmp00(i1:i2), *) char12
                 call util_unitstate(char12, inunit, instate)
                 do imp = inunit(1), inunit(2)
                    if(imp <= 0 .or. imp > iulen) then
                       error_message = 'Error: invalid value for particle id in specifying charged sites'
                       call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    inele%flag_ele(imp + lunit2mp(1, iu) - 1) = .false.
                 end do
              end if
              i1 = 0
              isw = 0
           end if
        end do
     else if(ctmp00(1:10) == 'CHARGE_ADD') then
        write (lunout, '(2a)') '---reading CHARGE_ADD: ', trim(ctmp00)
        isw = 0
        i1 = 0
        i2 = 0
        ifield = 0
        do k = 1, CARRAY_MXCOLM
           if(cwkinp(iline)(k:k) /= ' ' .and. &
                cwkinp(iline)(k:k) /= ',') then
              if(isw == 0) then
                 i1 = k
                 isw = 1
              end if
              i2 = k
           else if (isw == 1) then
              ifield = ifield + 1
              if(ifield == 2) then
                 read (ctmp00(i1:i2), *) char12
                 call util_unitstate(char12, inunit, instate)
                 iu = inunit(1)   !only consider real unit
                 if(iu <= 0 .or. iu > nunit_real) then
                    error_message = 'Error: invalid value for unit id in specifying charged sites'
                    call util_error(ERROR%STOP_ALL, error_message)
                 end if
                 iulen = lunit2mp(2, iu) - lunit2mp(1, iu) + 1
                 if(.not. icalc(iu)) then
                    write(error_message, '(a, i3, a)') 'Warning: unit ', iu,&
                       ' is not assigned electrostatic interaction '&
                       &//' but is specified the charged sites. Ignore'
                    call util_error(ERROR%WARN_ALL, error_message)
                    cycle loop_line 
                 end if
              else if(ifield == 3) then
                 read (ctmp00(i1:i2), *) char12
                 if(char12(1:1) /= 'u') then
                    error_message = 'Error: bad format for CHARGE_ADD line.  Forget the symbol u before particle id?'
                    call util_error(ERROR%STOP_ALL, error_message)
                 end if
              else if(ifield >= 4) then
                 read (ctmp00(i1:i2), *) char12
                 call util_unitstate(char12, inunit, instate)
                 do imp = inunit(1), inunit(2)
                    if(imp <= 0 .or. imp > iulen) then
                       error_message = 'Error: invalid value for particle id in specifying charged sites'
                       call util_error(ERROR%STOP_ALL, error_message)
                    end if
                    inele%flag_ele(imp + lunit2mp(1, iu) - 1) = .true.
                 end do
              end if
              i1 = 0
              isw = 0
           end if
        end do
     else if(ctmp00(1:13) == 'CHARGE_CHANGE') then
        i_charge_change = i_charge_change + 1
        if (i_charge_change > MXCHARGECHANGE) then
           error_message = 'Error: should increase MXCHARGECHANGE'
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        read (ctmp00, *) header, read_charge_change_imp, read_charge_change_value
        write (lunout, '(2a, i10, g10.3)') '---reading CHARGE_CHANGE: ', header, read_charge_change_imp, read_charge_change_value
        inele%charge_change_imp(i_charge_change) = read_charge_change_imp
        inele%charge_change_value(i_charge_change) = read_charge_change_value
     end if

  
  end do loop_line

  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     if (ctmp00(1:10) == 'DEL_CHARGE') then
        do icol = 12, CARRAY_MXCOLM
           if(ctmp00(icol:icol) == ')') exit
        end do

        read (ctmp00(12:icol-1), *) char12
        write (lunout, '(2a)') '---reading DEL_CHARGE: ', trim(char12)
        call util_unitstate(char12, inunit, instate)

        do imp = inunit(1), inunit(2)
           inele%flag_ele(imp) = .false.
        end do
     end if
  end do

  inele%n_charge_change = i_charge_change

  ! check
  if (inele%ionic_strength < ZERO_JUDGE) then
     error_message = 'Error: ionic_strength is too small'
     call util_error(ERROR%STOP_ALL, error_message)
  endif

  if(inele%i_diele == 0) then
     write (lunout, *) 'using constant value for dielectric constant'
     if (inele%diele_water > INVALID_JUDGE) then
        error_message = 'Error: invalid value for inele%diele_water'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  else if (inele%i_diele == 1) then
     write (lunout, *) "using dielectric constant as a function of temperature (Malmberg and Maryott, 1956)" 
  else if (inele%i_diele == 2) then
     write (lunout, *) "using dielectric constant as a function of temperature (used by Denesyuk)" 
  else
     error_message = 'Error: invalid value for inele%i_diele'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  if(inele%i_charge == 0) then
     continue  ! default
  else if (inele%i_charge == 1) then
     write (lunout, *) "using charge values based on Manning's theory"
  else
     error_message = 'Error: invalid value for inele%i_charge'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  if(inele%i_function_form == 0) then
     continue  ! default
  else if (inele%i_function_form == 1) then
     continue
     !write (lunout, *) "Coulomb potential is enabled."
  else if (inele%i_function_form == 2) then
     continue
     !write (lunout, *) "Coulomb potential is enabled. (Ewald method)"
  else if (inele%i_function_form == 3) then
     continue
     !write (lunout, *) "Coulomb potential is enabled. (Brute-force)"
  else
     error_message = 'Error: invalid value for inele%i_function_form'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  if (inele%i_function_form == 2) then ! Ewald summation
     if (inele%ewld_alpha > INVALID_JUDGE) then
        error_message = 'Error: invalid value for inele%ewld_alpha'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
     if (inele%ewld_hmax > INVALID_JUDGE) then
        error_message = 'Error: invalid value for inele%ewld_hmax'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
     if (inele%ewld_dipole < 0) then
        error_message = 'Error: invalid value for inele%ewld_dipole'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  endif

  if (inele%i_calc_method /= 0) then
     error_message = 'Error: invalid value for inele%i_calc_method (only 0 is valid in the current code)'
     call util_error(ERROR%STOP_ALL, error_message)
  endif

  ! output for .data
  !if (nlines /= 0) then
  if (any(inele%flag_ele)) then
     do aaid = 1, SEQID%CL
        write(lunout,'(a,1x,a,1x)',advance='NO') 'CHARGE_TYPE',cfunc_id2seq(aaid)
        write(lunout,'(f10.5)') inele%coef_charge_type(aaid)
     enddo
     write(lunout,'(a)') '>>>>'
  endif

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inele, inele%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inpmf, inpmf%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_electrostatic
