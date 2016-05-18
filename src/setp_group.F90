! setp_group
!> @brief Read the "<<<< group" field from .inp file.

subroutine setp_group()

  use const_maxsize
  use const_index
  use var_io,    only : infile, outfile
  use var_struct, only : nmp_all, cmass_mp, grp

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  integer :: icol, jcol, kcol
  integer :: luninp, lunout
  integer :: iline, nlines
  integer :: igrp
  integer :: h, iswp
  integer :: imp, mp1, mp2, nmp
  real(PREC) :: mass
  logical :: flg_swp
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  grp%ngrp = 0
  grp%nmp(:) = 0
  grp%implist(:,:) = 0
  grp%mass_fract(:,:) = 0.0e0_PREC

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'group           ', kfind, &
                     CARRAY_MXLINE, nlines, cwkinp)

  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     if(ctmp00(1:6) == 'GROUP(') then

        icol = 7        
        do while (ctmp00(icol:icol) /= ")")
           icol = icol+1
           if (icol == CARRAY_MXCOLM) then
              error_message = 'Error: format is wrong in line, ' // ctmp00
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        enddo

        read(ctmp00(7:icol-1),*) igrp
        if (igrp > MXGRP) then
           error_message = 'Error: Please increase MXGRP in const_maxsize.F90'
           call util_error(ERROR%STOP_ALL, error_message)
        else if (igrp > grp%ngrp) then
           grp%ngrp = igrp
        endif

        do while (ctmp00(icol:icol) /= "(")
           icol = icol+1
           if (icol == CARRAY_MXCOLM) then
              error_message = 'Error: format is wrong in line, ' // ctmp00
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        enddo

        nmp = 0
        icol = icol+1  ! Beginning
        kcol = 0       ! Position of "-" if appeared, otherwise 0

        ! jcol searches the position of "/" or ")"
        do jcol = icol+1, CARRAY_MXCOLM

           if (ctmp00(jcol:jcol) == "-") then
              kcol = jcol

           else if (ctmp00(jcol:jcol) == "/" .or. ctmp00(jcol:jcol) == ")") then

              ! /mp1/ style
              if (kcol == 0) then
                 read(ctmp00(icol:jcol-1), *) mp1
                 if (mp1 <= 0 .or. mp1 > nmp_all) then
                    error_message = &
                    'Error: format is wrong, or specified mp does not exist in line, ' // ctmp00
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif

                 nmp = nmp + 1
                 if (nmp > MXMPGRP) then
                    error_message = 'Error: Please increase MXMPGRP in const_maxsize.F90'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif
                 grp%implist(nmp, igrp) = mp1

              ! /mp1-mp2/ style
              else
                 read(ctmp00(icol:kcol-1),  *) mp1
                 read(ctmp00(kcol+1:jcol-1),*) mp2
                 if (mp1<=0 .or. mp2<=0 .or. mp1>nmp_all .or. mp2>nmp_all .or. mp1>mp2) then
                    error_message = &
                    'Error: format is wrong, or specified mp does not exist in line, ' // ctmp00
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif

                 do imp = mp1, mp2
                    nmp = nmp + 1
                    if (nmp > MXMPGRP) then
                       error_message = 'Error: Please increase MXMPGRP in const_maxsize.F90'
                       call util_error(ERROR%STOP_ALL, error_message)
                    endif
                    grp%implist(nmp, igrp) = imp
                 enddo
              endif
              
              ! Finish
              if (ctmp00(jcol:jcol) == ")") then
                 exit
              else if (jcol == CARRAY_MXCOLM) then
                 error_message = 'Error: format is wrong in line, ' // ctmp00
                 call util_error(ERROR%STOP_ALL, error_message)
              endif

              icol = jcol + 1
              kcol = 0

           else if (jcol == CARRAY_MXCOLM) then
              error_message = 'Error: format is wrong in line, ' // ctmp00
              call util_error(ERROR%STOP_ALL, error_message)
           end if
        end do

        grp%nmp(igrp) = nmp

        ! Sort (comb11)
        h = nmp
        flg_swp = .false.
        do while (h > 1 .or. flg_swp)
           h = (h * 10) / 13
           if (h < 1) then
              h = 1
           else if (h == 9 .or. h == 10) then 
              h = 11
           endif
           flg_swp = .false.
           do imp = 1, nmp-h
              if (grp%implist(imp,igrp) > grp%implist(imp+h,igrp)) then
                 iswp = grp%implist(imp,igrp)
                 grp%implist(imp,igrp) = grp%implist(imp+h,igrp)
                 grp%implist(imp+h,igrp) = iswp
                 flg_swp = .true.
              endif
           enddo
        enddo

        ! Check non-redundancy
        do imp = 1, nmp-1
           if (grp%implist(imp,igrp) == grp%implist(imp+1,igrp)) then
              error_message = 'Error: Same mp-ID specified more than twise in line, ' // ctmp00
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        enddo

     else
        !call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
        !do iequa = 1, nequat
        !   cvalue = 'i_lower_bound'
        !   call ukoto_ivalue2(lunout, csides(1, iequa), &
        !        inmisc%i_lower_bound, cvalue)
        !end do
     end if
  end do

  do igrp = 1, grp%ngrp
     write(lunout,'(a,i0,a,i0,a)',advance='NO') 'GROUP(',igrp,') #mp=',grp%nmp(igrp),' : '
     do imp = 1, grp%nmp(igrp)-1
        write(lunout,'(i0,a)',advance='NO') grp%implist(imp,igrp),' '
     enddo
     write(lunout,'(i0,a)') grp%implist(grp%nmp(igrp),igrp)
  enddo

#ifdef MPI_PAR
  end if
  call MPI_Bcast(grp, grp%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif

  !calc mass fraction (This is used in caluculation for center of mass.)
  do igrp = 1, grp%ngrp
     mass = 0.0e0_PREC
     do imp = 1, grp%nmp(igrp)
        mass = mass + cmass_mp(grp%implist(imp,igrp))
     enddo
     do imp = 1, grp%nmp(igrp)
        grp%mass_fract(imp,igrp) = cmass_mp(grp%implist(imp,igrp)) / mass
     enddo
  enddo

endsubroutine setp_group
