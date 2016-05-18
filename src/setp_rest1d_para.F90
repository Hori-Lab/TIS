! setp_rest1d_para
!> @brief This subroutine is to read and set the parameters for 1D-restraint simulation, &
!>        in which some mass-points are restrained to some positions on 1D-axis          &
!>        by harmonic springs.


! **********************************************************************
subroutine setp_rest1d_para()

  use const_maxsize
  use const_index
  use var_io, only : infile, outfile
  use var_setp, only : inmisc

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  ! local variables
  integer :: luninp, lunout
  integer :: iline, nlines
  integer :: irest, irest_center
  integer :: icol, jcol, imp
  integer :: inunit(2), instate
  character(4) :: kfind
  character(6) :: char6
  character(13) :: char12
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data
  inmisc%nrest1d = 0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'rest1d_para     ', kfind,&
                     CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "rest1d_para" field in set_rest1d_para'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  irest = 0
  irest_center = 0
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)

     if (ctmp00(1:13) == 'REST1D_CENTER') then
        irest_center = irest_center + 1
        inmisc%nrest1d_center_mp(irest_center) = 0
        icol = 15        
        do jcol = 8, CARRAY_MXCOLM
           if (ctmp00(jcol:jcol) == "/" .or. ctmp00(jcol:jcol) == ")") then
              read(ctmp00(icol:jcol-1), *) char12
              call util_unitstate(char12, inunit, instate)
              do imp = inunit(1), inunit(2)
                 inmisc%nrest1d_center_mp(irest_center) = inmisc%nrest1d_center_mp(irest_center) + 1
                 inmisc%irest1d_center2mp( irest_center, inmisc%nrest1d_center_mp(irest_center) ) = imp
              end do
              icol = jcol + 1
              if (ctmp00(jcol:jcol) == ")") then
                 exit;
              end if
           end if
        end do
        read(ctmp00(icol:CARRAY_MXCOLM), *) inmisc%coef_rest1d_center(irest_center),   inmisc%irest1d_center_mp_sa(irest_center), &
                                            inmisc%irest1d_center_mp_sb(irest_center), inmisc%rest1d_center_s0(irest_center)
        inmisc%rest1d_center_init_flag(irest_center) = 1
     else if(ctmp00(1:6) == 'REST1D') then
        irest = irest + 1
        read(ctmp00, *) char6, inmisc%irest1d2mp(irest), inmisc%coef_rest1d(irest),&
                               inmisc%irest1d_mp_sa(irest), inmisc%irest1d_mp_sb(irest), &
                               inmisc%rest1d_s0(irest)
     end if
  end do
  inmisc%nrest1d = irest
  inmisc%nrest1d_center = irest_center

  if(inmisc%nrest1d > MXREST1D) then
     error_message = 'Error: should increase MXREST1D'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  write(lunout,*) '<<<< rest1d_para'
  do irest = 1, inmisc%nrest1d
     write(lunout,*) 'REST1D:',irest
     write(lunout,*) 'mp_target=',inmisc%irest1d2mp(irest)
     write(lunout,*) 'coef=',inmisc%coef_rest1d(irest)
     write(lunout,*) 'mp_sa=',inmisc%irest1d_mp_sa(irest)
     write(lunout,*) 'mp_sb=',inmisc%irest1d_mp_sb(irest)
     write(lunout,*) 's0=',inmisc%rest1d_s0(irest)
  end do

  do irest_center = 1, inmisc%nrest1d_center
     write(lunout,*) 'REST1D_CENTER:', irest_center
     write(lunout,*) 'mp_target='
     do imp = 1, inmisc%nrest1d_center_mp(irest_center)
        write(lunout, *) inmisc%irest1d_center2mp(irest_center, imp)
     end do
     write(lunout,*) 'coef=',inmisc%coef_rest1d_center(irest_center)
     write(lunout,*) 'mp_sa=',inmisc%irest1d_center_mp_sa(irest_center)
     write(lunout,*) 'mp_sb=',inmisc%irest1d_center_mp_sb(irest_center)
     write(lunout,*) 's0=',inmisc%rest1d_center_s0(irest_center)
  end do
   
  write(lunout,*) '>>>>'

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inmisc, inmisc%sz ,MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_rest1d_para
