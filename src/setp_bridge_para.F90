! setp_bridge_para
!> @brief Read parameters for bridge option

! **********************************************************************
subroutine setp_bridge_para()

  use const_maxsize
  use const_index
  use var_inp, only : infile, outfile, i_run_mode
  use var_setp, only : inmisc

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  ! intent(out) :: 

  ! --------------------------------------------------------------------
  ! local variables
  integer :: ioerr
  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  integer :: ibrid, ibrid_com, ibrid_ppr
  real(PREC) :: rzero
  character(6) :: char6
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data
  inmisc%nbrid = 0
  inmisc%nbrid_com = 0
  inmisc%nbrid_ppr = 0
  inmisc%i_lower_bound = 0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'bridge_para     ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "bridge_para" field in setp_bridge_para'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  ibrid = 0
  ibrid_com = 0
  ibrid_ppr = 0
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     if (ctmp00(1:13) == 'BRIDGE_CENTER') then
        ibrid_com = ibrid_com+1
        if(ibrid_com > MXBRIDGE) then
           error_message = 'Error: should increase MXBRIDGE'
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        read (ctmp00, *) char6, &
             inmisc%ibrid_com2grp(1,ibrid_com), inmisc%ibrid_com2grp(2,ibrid_com), &
             inmisc%coef_brid_com(ibrid_com), inmisc%brid_com_dist(ibrid_com)
        write (lunout, '(1a, 2i10, 2g10.3)') &
             '---reading bridge residues: BRIDGE_CENTER', &
             inmisc%ibrid_com2grp(1,ibrid_com), inmisc%ibrid_com2grp(2,ibrid_com), &
             inmisc%coef_brid_com(ibrid_com), inmisc%brid_com_dist(ibrid_com)

     else if(ctmp00(1:10) == 'BRIDGE_PPR') then
        if (i_run_mode == RUN%REPLICA) then
           error_message = 'Error: BRIDGE_PPR does not work in REM simulation.'
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        ibrid_ppr = ibrid_ppr+1
        ibrid_com = ibrid_com+1
        inmisc%ibrid_ppr_com(ibrid_ppr) = ibrid_com
        inmisc%ibrid_ppr_cycl(ibrid_ppr) = 0  ! to check whether reading is OK
        inmisc%ibrid_ppr_opt(ibrid_ppr) = 0   ! default: no option
        if(ibrid_ppr > MXBRIDGE .or. ibrid_com > MXBRIDGE) then
           error_message = 'Error: should increase MXBRIDGE'
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        read (ctmp00, *,iostat=ioerr) char6, &
             inmisc%ibrid_com2grp(1,ibrid_com), inmisc%ibrid_com2grp(2,ibrid_com), &
             inmisc%ibrid_ppr_gid_r(1,ibrid_ppr), inmisc%ibrid_ppr_gid_r(2,ibrid_ppr), &
             inmisc%coef_brid_com(ibrid_com),&
             inmisc%brid_ppr_rmin(ibrid_ppr), inmisc%brid_ppr_rcut(ibrid_ppr), &
             inmisc%brid_ppr_rmax(ibrid_ppr), rzero,&
             inmisc%ibrid_ppr_cycl(ibrid_ppr), inmisc%ibrid_ppr_opt(ibrid_ppr)
        ! ioerr can be 0 because ibrid_ppr_opt is optional.
        ! Thus ibrid_ppr_cycl shoule be checked.
        if (ioerr > 0 .or. inmisc%ibrid_ppr_cycl(ibrid_ppr) <= 0) then
           error_message = 'Error: Bad format in BRIDGE_PPR line.'
           call util_error(ERROR%STOP_ALL, error_message)
        endif

        write (lunout, '(1a, 4i4, 5g8.3, i11, i2)') &
             '---reading bridge residues: BRIDGE_PPR', &
             inmisc%ibrid_com2grp(1,ibrid_com), inmisc%ibrid_com2grp(2,ibrid_com), &
             inmisc%ibrid_ppr_gid_r(1,ibrid_com), inmisc%ibrid_ppr_gid_r(2,ibrid_com), &
             inmisc%coef_brid_com(ibrid_com),&
             inmisc%brid_ppr_rmin(ibrid_ppr), inmisc%brid_ppr_rcut(ibrid_ppr), &
             inmisc%brid_ppr_rmax(ibrid_ppr), rzero,&
             inmisc%ibrid_ppr_cycl(ibrid_ppr), inmisc%ibrid_ppr_opt(ibrid_ppr)

        inmisc%brid_ppr_rzero2 = rzero ** 2

     else if(ctmp00(1:6) == 'BRIDGE') then
        ibrid = ibrid + 1
        if(ibrid > MXBRIDGE) then
           error_message = 'Error: should increase MXBRIDGE'
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        read (ctmp00, *) char6, inmisc%ibrid2mp(1, ibrid), &
             inmisc%ibrid2mp(2, ibrid), &
             inmisc%coef_brid(ibrid), inmisc%brid_dist(ibrid)
        write (lunout, '(2a, 2i10, 2g10.3)') '---reading bridge residues: ', &
             char6, inmisc%ibrid2mp(1, ibrid), &
             inmisc%ibrid2mp(2, ibrid), &
             inmisc%coef_brid(ibrid), inmisc%brid_dist(ibrid)
     else
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
        
        do iequa = 1, nequat
           cvalue = 'i_lower_bound'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inmisc%i_lower_bound, cvalue)
        end do

     end if
  end do
  inmisc%nbrid = ibrid 
  inmisc%nbrid_com = ibrid_com
  inmisc%nbrid_ppr = ibrid_ppr
  
#ifdef MPI_PAR
  end if

  call MPI_Bcast(inmisc, inmisc%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_bridge_para
