! mloop_ninfo_inp
!> @brief Read the filename of ".ninfo"; called by mloop_nativeinfo

subroutine mloop_ninfo_inp(istep_sim, i_ninfo_type, inat_unit, cnat_fname)
  
  use const_maxsize
  use const_index
  use var_inp, only : infile, outfile, ius2unit

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  integer, intent(in) :: istep_sim
  integer, intent(out) :: i_ninfo_type
  integer, intent(out) :: inat_unit(MXUNIT, MXUNIT)
  character(CARRAY_MXFILE), intent(out) :: cnat_fname(MXUNIT*MXUNIT)

  ! --------------------------------------------------------------------
  ! local variables
  integer :: iu, ju, npath
  integer :: inat, icol, jcol, iunit, junit, num_file(2)
  integer :: inunit(2), jnunit(2), instate, instate2
  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat

  character(CARRAY_MXFILE) :: path_natinfo
  character(CARRAY_MXFILE) :: cnat
  character(12) :: char12
  character(16) :: char16
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: char00
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MXCOLM) :: ctmp02
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  iunit = 0
  inat_unit(1:MXUNIT, 1:MXUNIT) = 0
  do inat = 1, MXUNIT*MXUNIT
     cnat_fname(inat) = ""
  end do

  ! -------------------------------------------------------
  ! read native information from luninp 
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  path_natinfo = "./ninfo"
  rewind(luninp)
  call ukoto_uiread2(luninp, 6, 'filenames       ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  do iline = 1, nlines
     call ukoto_uiequa2(6, cwkinp(iline), nequat, csides)
     do iequa = 1, nequat
        if(csides(1, iequa) == 'path_natinfo') then
           path_natinfo = csides(2, iequa)
        end if
     end do
  end do
  npath = index(path_natinfo, ' ')

  rewind(luninp)
  if((istep_sim / 10) >= 1) then 
     write(char16, '(''native_info_sim'',i2)') istep_sim
  else
     write(char16, '(''native_info_sim'',i1)') istep_sim
  end if
  call ukoto_uiread2(luninp, lunout, char16, kfind, &
       CARRAY_MXLINE, nlines, cwkinp)                                 
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "native_info_sim" field in mloop_ninfo_inp'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  num_file(1:2) = 0
  do iline = 1, nlines
     call ukoto_uiequa2(6, cwkinp(iline), nequat, csides)
     do iequa = 1, nequat
        if(csides(1, iequa) == 'path') then
           path_natinfo = csides(2, iequa)
           npath = index(path_natinfo, ' ')
        end if
     end do
  end do

  i_ninfo_type = NATIVEINFO%ONE_BY_ONE
  do iline = 1, nlines 
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     char00 = cwkinp(iline)
     
     if(char00(1:14) == 'NINFO(all/all)') then
        inat = 1
        i_ninfo_type = NATIVEINFO%ALL_IN_ONE

        do icol = 15, CARRAY_MXCOLM
           if(char00(icol:icol) /= ' ') exit
        end do
        cnat = char00(icol:CARRAY_MXCOLM)
        if(path_natinfo /= '') then
           cnat_fname(inat) = path_natinfo(1:npath-1)//'/'//cnat
        else
           cnat_fname(inat) = cnat
        end if
        num_file(1) = num_file(1) + 1
     end if
  end do

  do iline = 1, nlines 
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     char00 = cwkinp(iline)
     
     if(char00(1:5) == 'NINFO' .and. char00(1:14) /= 'NINFO(all/all)') then

        do icol = 7, CARRAY_MXCOLM
           if(char00(icol:icol) == '/') exit
        end do
        do jcol = 7, CARRAY_MXCOLM
           if(char00(jcol:jcol) == ')') exit
        end do
        
        read (char00(7:icol-1), '(a)') char12
        call util_unitstate(char12, inunit, instate)
        read (char00(icol+1:jcol-1), '(a)') char12
        call util_unitstate(char12, jnunit, instate2)
        
        read (char00(jcol+1:CARRAY_MXCOLM), *) inat
        write (lunout, '(3a, i4)') '---reading nativeinfo: ', char00(1:jcol), ' ', inat
        do iu = inunit(1), inunit(2)
           do ju = jnunit(1), jnunit(2)
              iunit = ius2unit(iu, instate)
              junit = ius2unit(ju, instate2)
              inat_unit(iunit, junit) = inat
              iunit = iunit + 1 
           end do
        end do
     end if
     
     do iequa = 1, nequat
        ctmp02 = csides(1, iequa)
        if(ctmp02(1:1) >= '1' .and. ctmp02(1:1) <= '9') then
           read (csides(1, iequa), *) inat
           cnat = csides(2, iequa)
           write (lunout, '(a, i4, 2a)') '---reading nativeinfo file name: ', inat, ' ', trim(cnat)
           if(path_natinfo /= '') then
              cnat_fname(inat) = path_natinfo(1:npath-1)//'/'//cnat
           else
              cnat_fname(inat) = cnat
           end if
           num_file(2) = num_file(2) + 1
        end if
     end do
  end do
  
  if(num_file(1) >= 1 .and. num_file(2) >= 1) then
     error_message = 'Error: reading native_info file type should be 1 in "native_info_sim" field'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast (i_ninfo_type,        1,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (inat_unit,MXUNIT*MXUNIT,MPI_INTEGER,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (cnat_fname,MXUNIT*MXUNIT*CARRAY_MXFILE,MPI_CHARACTER ,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine mloop_ninfo_inp
