!setp_pulling_para
!> @brief Reads "pulling_para" field and stores the data into  &
!>        "inmisc" struct.

subroutine setp_pulling_para()

  use const_maxsize
  use const_index
  use const_physical
  use var_inp, only : infile, outfile
  use var_setp, only : inmisc
  use var_struct, only: nmp_all, lunit2mp
  use var_replica, only: flg_rep, n_replica_all, rep2val, inrep, lab2val

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  ! intent(out) :: 

  ! --------------------------------------------------------------------
  ! local variables
  integer :: idx, imp, icol, iunit, grep
  integer :: luninp, lunout
  integer :: iline, nlines, iequat, nequat
  integer :: ipull, iunravel
  integer :: inunit(2), instate
  real(PREC) :: factor, pullforce
  real(PREC) :: local_pull_xyz(3)

  character(7)  :: char7
  character(11) :: char11
  character(12) :: char12
  character(4)  :: kfind
  character(CARRAY_MXCOLM)  :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MXCOLM)  :: cvalue
  character(CARRAY_MXCOLM)  :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: cdummy
  character(CARRAY_MXFILE) :: filename
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data
  inmisc%npull = 0
  inmisc%coef_pull(1:MXPULLING) = -1.0
  inmisc%pu_xyz(1:3, 1:MXPULLING) = 0.0

  inmisc%ipull_force_unit = 0  ! optional, (default:0 kcal/mol/A)

  if (flg_rep(REPTYPE%PULL)) then
     !! When REM is enabled, npull_unravel is already set in inp_replica_para
     do grep = 1, n_replica_all
        pullforce = rep2val(grep, REPTYPE%PULL)
        inmisc%pull_unravel_xyz(:, 1:inrep%n_pull, grep) = pullforce * inrep%pull_direction(:, 1:inrep%n_pull)
     enddo
     iunravel = inmisc%npull_unravel
  else
     iunravel = 0
  endif
  inmisc%ipull_unravel2mp(1:2, iunravel+1:MXPULLING) = 0
  inmisc%pull_unravel_xyz(1:3, iunravel+1:MXPULLING, 1:MXREPLICA) = 0.0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'pulling_para    ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "pulling_para" field in set_pulling_para'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  ipull = 0
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     if (ctmp00(1:12) == 'PULL_CF_FILE') then
        idx = index(trim(ctmp00), " ", back=.true.) 
        read (ctmp00(1:idx), *) char12, factor
        write (lunout, '(a, g10.3)') '---reading PULL_CF_FILE: ', factor
        filename = trim(ctmp00(idx+1:))
        call read_pull(ipull, factor, filename)

     else if(ctmp00(1:7) == 'PULL_CF') then
        ipull = ipull + 1
        read (ctmp00, *) char7, inmisc%ipull2mp(ipull), &
             inmisc%pull_xyz(1, ipull), inmisc%pull_xyz(2, ipull), &
             inmisc%pull_xyz(3, ipull)
        write (lunout, '(a, i10, 3g10.3)') '---reading PULL_CF: ', inmisc%ipull2mp(ipull), &
             inmisc%pull_xyz(1, ipull), inmisc%pull_xyz(2, ipull), &
             inmisc%pull_xyz(3, ipull)

     else if(ctmp00(1:10) == 'PULL_MP_CF') then
        
        do icol = 12, CARRAY_MXCOLM
           if(ctmp00(icol:icol) == ')') exit
        end do

        read (ctmp00(12:icol-1), *) char12
        write (lunout, '(2a)') '---reading PULL_MP_CF: ', trim(char12)
        call util_unitstate(char12, inunit, instate)

        read (ctmp00(icol+1:), *) local_pull_xyz(1), local_pull_xyz(2), local_pull_xyz(3)
        
        do imp = inunit(1), inunit(2)
           ipull = ipull + 1
           inmisc%ipull2mp(ipull) = imp
           inmisc%pull_xyz(1, ipull) = local_pull_xyz(1)
           inmisc%pull_xyz(2, ipull) = local_pull_xyz(2)
           inmisc%pull_xyz(3, ipull) = local_pull_xyz(3)
        end do
        
     else if(ctmp00(1:12) == 'PULL_UNIT_CF') then
        
        do icol = 14, CARRAY_MXCOLM
           if(ctmp00(icol:icol) == ')') exit
        end do

        read (ctmp00(14:icol-1), *) char12
        write (lunout, '(2a)') '---reading PULL_UNIT_CF: ', trim(char12)
        call util_unitstate(char12, inunit, instate)

        read (ctmp00(icol+1:), *) local_pull_xyz(1), local_pull_xyz(2), local_pull_xyz(3)

        do iunit = inunit(1), inunit(2)
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
              ipull = ipull + 1
              inmisc%ipull2mp(ipull) = imp
              inmisc%pull_xyz(1, ipull) = local_pull_xyz(1)
              inmisc%pull_xyz(2, ipull) = local_pull_xyz(2)
              inmisc%pull_xyz(3, ipull) = local_pull_xyz(3)
           end do
        end do

     else if(ctmp00(1:11) == 'PULL_ALL_CF') then
        read (ctmp00, *) char11, &
             local_pull_xyz(1), local_pull_xyz(2), local_pull_xyz(3)
        write (lunout, '(a, 3g10.3)') '---reading PULL_ALL_CF: ', &
             local_pull_xyz(1), local_pull_xyz(2), local_pull_xyz(3)

        do imp = 1, nmp_all
           ipull = ipull + 1
           inmisc%ipull2mp(ipull) = imp
           inmisc%pull_xyz(1, ipull) = local_pull_xyz(1)
           inmisc%pull_xyz(2, ipull) = local_pull_xyz(2)
           inmisc%pull_xyz(3, ipull) = local_pull_xyz(3)
        end do

     else if(ctmp00(1:7) == 'PULL_CV') then
        ipull = ipull + 1
        read (ctmp00, *) char7, inmisc%ipull2mp(ipull), &
             inmisc%coef_pull(ipull), &
             inmisc%pull_xyz(1, ipull), inmisc%pull_xyz(2, ipull), &
             inmisc%pull_xyz(3, ipull), &
             inmisc%pu_xyz(1, ipull), inmisc%pu_xyz(2, ipull), &
             inmisc%pu_xyz(3, ipull)
        write (lunout, '(a, i10, 7g10.3)') '---reading PULL_CV: ', inmisc%ipull2mp(ipull), &
             inmisc%coef_pull(ipull), &
             inmisc%pull_xyz(1, ipull), inmisc%pull_xyz(2, ipull), &
             inmisc%pull_xyz(3, ipull), &
             inmisc%pu_xyz(1, ipull), inmisc%pu_xyz(2, ipull), &
             inmisc%pu_xyz(3, ipull)

     else if(ctmp00(1:15) == 'PULL_UNRAVEL_CF') then
        iunravel = iunravel + 1
        read (ctmp00, *) cdummy, inmisc%ipull_unravel2mp(1,iunravel), &
                                 inmisc%ipull_unravel2mp(2,iunravel), &
                                 inmisc%pull_unravel_xyz(1,iunravel,1), &
                                 inmisc%pull_unravel_xyz(2,iunravel,1), &
                                 inmisc%pull_unravel_xyz(3,iunravel,1)
        write (lunout, '(a, i5, 1x, i5, 3g10.3)') '---reading PULL_CV: ', &
                                 inmisc%ipull_unravel2mp(1,iunravel), &
                                 inmisc%ipull_unravel2mp(2,iunravel), &
                                 inmisc%pull_unravel_xyz(1,iunravel,1), &
                                 inmisc%pull_unravel_xyz(2,iunravel,1), &
                                 inmisc%pull_unravel_xyz(3,iunravel,1)
     else
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
          
        do iequat = 1, nequat
           cvalue = 'i_force_unit'
           call ukoto_ivalue2(lunout, csides(1, iequat), &
                              inmisc%ipull_force_unit, cvalue)
        enddo
     end if
  end do
  inmisc%npull = ipull
  inmisc%npull_unravel = iunravel
  
  if(inmisc%npull > MXPULLING) then
     error_message = 'Error: should increase MXPULLING'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  if(inmisc%npull_unravel > MXPULLING) then
     error_message = 'Error: should increase MXPULLING'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! output to datafile
  write(lunout,'(a)') '<<<< pulling_para'
  do ipull = 1, inmisc%npull
     if (inmisc%coef_pull(ipull) <= 0.0e0_PREC) then
        write(lunout, '(a)', advance='no') 'PULL_CF '
        write(lunout, '(i6,3(1x,f11.6))')  inmisc%ipull2mp(ipull), inmisc%pull_xyz(1:3, ipull)
     else
        write(lunout, '(a)', advance='no') 'PULL_CV '
        write(lunout, '(i6,7(1x,f11.6))') inmisc%ipull2mp(ipull), inmisc%coef_pull(ipull), &
                                        inmisc%pull_xyz(1:3, ipull), inmisc%pu_xyz(1:3, ipull)
     endif
  enddo
  do iunravel = 1, inmisc%npull_unravel
     write(lunout, '(a)', advance='no') 'PULL_UNRAVEL '
     write(lunout, '(i6,1x,i6,3(1x,f11.6))') inmisc%ipull_unravel2mp(1,iunravel), &
                                           inmisc%ipull_unravel2mp(2,iunravel), &
                                           inmisc%pull_unravel_xyz(1:3,iunravel,1)
  enddo
  write(lunout,'(a)') '>>>>'

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inmisc, inmisc%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif

  ! When force values are given in pN, convert them to kcal/mol/A
  if (inmisc%ipull_force_unit == 1) then
     do ipull = 1, inmisc%npull
        if (inmisc%coef_pull(ipull) < 0.0) then
           inmisc%pull_xyz(1:3, ipull) = inmisc%pull_xyz(1:3, ipull) &
                                        * N_AVO * 1.0e-12 * 1.0e-10 / JOUL2KCAL
        endif
     enddo

     inmisc%pull_unravel_xyz(:,:,:) = inmisc%pull_unravel_xyz(:,:,:) &
                                     * N_AVO * 1.0e-12 * 1.0e-10 / JOUL2KCAL
     if (flg_rep(REPTYPE%PULL)) then
        lab2val(1:n_replica_all, REPTYPE%PULL) = lab2val(1:n_replica_all, REPTYPE%PULL) &
                                                * N_AVO * 1.0e-12 * 1.0e-10 / JOUL2KCAL
     endif
  endif

end subroutine setp_pulling_para
