!setp_hydrophobic
!> @brief Reads "<<<< hydrophobic" field to redefine parameters
!>        for hydrophobic interaction; this subroutine also
!>        sets the hydrophobic sites

subroutine setp_hydrophobic()

  use const_maxsize
  use const_index
  use var_io, only : infile, outfile
  use var_struct, only : nunit_real, lunit2mp
  use var_setp, only : inhp, inmisc

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  ! intent(out) :: xbox, ybox, zbox
  
  ! --------------------------------------------------------------------
  ! local variables
  integer :: i, j, icol, jcol, isw, i1, i2, ifield, iunit, junit
  integer :: k, imp, aaid
  integer :: inunit(2), jnunit(2), instate, iu, iulen
  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  integer :: ncoor
  real(PREC) :: ncoormax
  real(PREC) :: x, chp 
  logical :: icalc(MXUNIT)

  character(3) :: char3
  character(12) :: char12
  character(50) :: char50
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MXCOLM) :: ctmp02
  character(CARRAY_MSG_ERROR) :: error_message
  
  integer :: ifunc_seq2id

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data
  inhp%flag_hp(1:MXMP) = .false.
  icalc(1:MXUNIT) = .false. 
  do iunit = 1, nunit_real
     do junit = iunit, nunit_real
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%HP)) then
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
              inhp%flag_hp(imp) = .true.
           end do
           do imp = lunit2mp(1, junit), lunit2mp(2, junit)
              inhp%flag_hp(imp) = .true.
           end do
           icalc(iunit) = .true.
           icalc(junit) = .true.
        end if
     end do
  end do

#ifdef MPI_PAR
  if (myrank == 0) then
#endif
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'hydrophobic     ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind == 'FIND') then
        write(lunout, '(71(1H*))')
        write(lunout, '(a)') '! following are parameters of hydrophobic interaction to be changed'
        do iline = 1, nlines
           ctmp00 = cwkinp(iline)
           if(ctmp00(1:8) == 'HPE_TYPE') then
              read (ctmp00, *) char50, char3, chp, ncoor, ncoormax
              write (lunout, '(3a, g10.3, i10, g10.3)') '---HPE_TYPE: ', trim(char50), char3, chp, ncoor, ncoormax
              if (char3 == 'OTH') then
                  aaid = 21
              else
                 aaid = ifunc_seq2id(char3)
              end if
              inhp%coefaa_para_hp(aaid) = chp
              inhp%ncoor_para_hp(aaid) = ncoor
              inhp%ncoormax_para_hp(aaid) = ncoormax
              write (lunout, '(a, a3, a, f10.5)') 'coef_aa_hp(', char3, ') = ', chp
              write (lunout, '(a, a3, a, i5)') 'rncoor_aa_hp(', char3, ') = ', ncoor
              write (lunout, '(a, a3, a, f10.5)') 'rncoormax_aa_hp(', char3, ') = ', ncoormax 
              cycle 
           end if
 
           call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
           do iequa = 1, nequat
              ctmp02 = csides(1, iequa)

              ! ---------------------------------------------------------------
              ! hydrophobic.para
              cvalue = 'coef_hp'
              if(csides(1, iequa) == cvalue) then
                 call ukoto_rvalue2(lunout, csides(1, iequa), &
                      inhp%coef_hp, cvalue)
                 write (lunout, '(a, f10.5)') 'coef_hp = ', inhp%coef_hp
              end if

              cvalue = 'coef_rho_hp'
              if(csides(1, iequa) == cvalue) then
                 call ukoto_rvalue2(lunout, csides(1, iequa), &
                      inhp%coef_rho_hp, cvalue)
                 write (lunout, '(a, f10.5)') 'coef_rho_hp = ', inhp%coef_rho_hp
              end if

              cvalue = 'rho_min_hp'
              if(csides(1, iequa) == cvalue) then
                 call ukoto_rvalue2(lunout, csides(1, iequa), &
                      inhp%rho_min_hp, cvalue)
                 write (lunout, '(a, f10.5)') 'rho_min_hp = ', inhp%rho_min_hp
              end if

              cvalue = 'coef_aa_hp'
              if(ctmp02(1:10) == cvalue) then
                 do icol = 1, CARRAY_MXCOLM
                    if(ctmp02(icol:icol) == ')') exit
                 end do
                 read (ctmp02(12:icol-1), *) char12
                 call util_unitstate(char12, inunit, instate)
                 cvalue = ctmp02(1:icol)
                 call ukoto_rvalue2(lunout, csides(1, iequa), &
                      x, cvalue)
                 do i = inunit(1), inunit(2)
                    inhp%coefaa_para_hp(i) = x
                 end do
                 write (lunout, '(a, a, a, f10.5)') 'Advanced user changes: ', ctmp02(1:icol), '=', x
              end if

              cvalue = 'rncoor_aa_hp'
              if(ctmp02(1:12) == cvalue) then
                 do icol = 1, CARRAY_MXCOLM
                    if(ctmp02(icol:icol) == ')') exit
                 end do
                 read (ctmp02(14:icol-1), *) char12
                 call util_unitstate(char12, inunit, instate)
                 cvalue = ctmp02(1:icol)
                 call ukoto_ivalue2(lunout, csides(1, iequa), &
                      i, cvalue)
                 do j = inunit(1), inunit(2)
                    inhp%ncoor_para_hp(j) = i
                 end do
                 write (lunout, '(a, a, a, i5)') 'Advanced user changes: ', ctmp02(1:icol), '=', i
              end if

              cvalue = 'rncoormax_aa_hp'
              if(ctmp02(1:15) == cvalue) then
                 do icol = 1, CARRAY_MXCOLM
                    if(ctmp02(icol:icol) == ')') exit
                 end do
                 read (ctmp02(17:icol-1), *) char12
                 call util_unitstate(char12, inunit, instate)
                 cvalue = ctmp02(1:icol)
                 call ukoto_rvalue2(lunout, csides(1, iequa), &
                      x, cvalue)
                 do j = inunit(1), inunit(2)
                    inhp%ncoormax_para_hp(j) = x
                 end do
                 write (lunout, '(a, a, a, f10.5)') 'Advanced user changes: ', ctmp02(1:icol), '=', x
              end if

              cvalue = 'cutoff_dmin_hp'
              if(ctmp02(1:14) == cvalue) then
                 do icol = 16, CARRAY_MXCOLM
                    if(ctmp02(icol:icol) == '/') exit
                 end do
                 do jcol = 16, CARRAY_MXCOLM
                    if(ctmp02(jcol:jcol) == ')') exit
                 end do
            
                 read (ctmp02(16:icol-1), *) char12
                 call util_unitstate(char12, inunit, instate)
                 read (ctmp02(icol+1:jcol-1), *) char12
                 call util_unitstate(char12, jnunit, instate)
                 
                 cvalue = ctmp02(1:jcol)
                 call ukoto_rvalue2(lunout, csides(1, iequa), &
                      x, cvalue)
                 do i = inunit(1), inunit(2)
                    do j = jnunit(1), jnunit(2)
                       inhp%cutoffdmin_para_hp(i,j) = x
                    end do
                 end do
                 write (lunout, '(a, a, a, f10.5)') 'Advanced user changes: ', ctmp02(1:jcol), '=', x
              end if

              cvalue = 'cutoff_dmax_hp'
              if(ctmp02(1:14) == cvalue) then
                 do icol = 16, CARRAY_MXCOLM
                    if(ctmp02(icol:icol) == '/') exit
                 end do
                 do jcol = 16, CARRAY_MXCOLM
                    if(ctmp02(jcol:jcol) == ')') exit
                 end do
            
                 read (ctmp02(16:icol-1), *) char12
                 call util_unitstate(char12, inunit, instate)
                 read (ctmp02(icol+1:jcol-1), *) char12
                 call util_unitstate(char12, jnunit, instate)
            
                 cvalue = ctmp02(1:jcol)
                 call ukoto_rvalue2(lunout, csides(1, iequa), &
                      x, cvalue)
                 do i = inunit(1), inunit(2)
                    do j = jnunit(1), jnunit(2)
                       inhp%cutoffdmax_para_hp(i,j) = x
                    end do
                 end do
                 write (lunout, '(a, a, a, f10.5)') 'Advanced user changes: ', ctmp02(1:jcol), '=', x
              end if
           end do
        end do
     end if  ! kfind

     ! following to read specific hydrophobic site
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'hydrophobic     ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind == 'FIND') then
       loop_line: do iline = 1, nlines
           ctmp00 = cwkinp(iline)
           if(ctmp00(1:7) == 'HPE_DEL') then
              write (lunout, '(a)') '---reading HPE_DEL: ', ctmp00
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
                           error_message = 'Error: invalid value for unit id in specifying hydrophobic sites'
                           call util_error(ERROR%STOP_ALL, error_message)
                        end if
                        iulen = lunit2mp(2, iu) - lunit2mp(1, iu) + 1
                        if(.not. icalc(iu)) then
                           write(error_message, '(a, i3, a)') 'Warning: unit ', iu, ' is not assigned hydrophobic interaction' &
                                   //' but is specified the hydrophobic sites. Ignore'
                           call util_error(ERROR%WARN_ALL, error_message)
                           cycle loop_line 
                        end if
                    else if(ifield == 3) then
                        read (ctmp00(i1:i2), *) char12
                        if(char12(1:1) /= 'u') then
                           error_message = 'Error: bad format for HPE_DEL line.  Forget the symbol u before particle id?'
                           call util_error(ERROR%STOP_ALL, error_message)
                        end if
                    else if(ifield >= 4) then
                        read (ctmp00(i1:i2), *) char12
                        call util_unitstate(char12, inunit, instate)
                        do imp = inunit(1), inunit(2)
                           if(imp <= 0 .or. imp > iulen) then
                              error_message = 'Error: invalid value for particle id in specifying hydrophobic sites'
                              call util_error(ERROR%STOP_ALL, error_message)
                           end if
                           inhp%flag_hp(imp + lunit2mp(1, iu) - 1) = .false.
                        end do
                    end if
                    i1 = 0
                    isw = 0
                 end if
              end do
           else if(ctmp00(1:7) == 'HPE_ADD') then
              write (lunout, '(2a)') '---reading HPE_ADD: ', ctmp00
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
                           error_message = 'Error: invalid value for unit id in specifying hydrophobic sites'
                           call util_error(ERROR%STOP_ALL, error_message)
                        end if
                        iulen = lunit2mp(2, iu) - lunit2mp(1, iu) + 1
                        if(.not. icalc(iu)) then
                           write(error_message, '(a, i3, a)') 'Warning: unit ', iu, &
                                     ' is not assigned hydrophobic interaction ' &
                                   //'but is specified the hydrophobic sites. Ignore'
                           call util_error(ERROR%WARN_ALL, error_message)
                           cycle loop_line 
                        end if
                    else if(ifield == 3) then
                        read (ctmp00(i1:i2), *) char12
                        if(char12(1:1) /= 'u') then
                           error_message = 'Error: bad format for HPE_ADD line.  Forget the symbol u before particle id?'
                           call util_error(ERROR%STOP_ALL, error_message)
                        end if
                    else if(ifield >= 4) then
                        read (ctmp00(i1:i2), *) char12
                        call util_unitstate(char12, inunit, instate)
                        do imp = inunit(1), inunit(2)
                           if(imp <= 0 .or. imp > iulen) then
                              error_message = 'Error: invalid value for particle id in specifying hydrophobic sites'
                              call util_error(ERROR%STOP_ALL, error_message)
                           end if
                           inhp%flag_hp(imp + lunit2mp(1, iu) - 1) = .true.
                        end do
                    end if
                    i1 = 0
                    isw = 0
                 end if
              end do
           end if
        end do loop_line
     end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast (inhp, inhp%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)

#endif

end subroutine setp_hydrophobic
