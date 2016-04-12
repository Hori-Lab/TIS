! setp_para_mgo
!> @brief 

! ***********************************************************************
subroutine setp_para_mgo()

  use const_maxsize
  use const_index
  use const_physical
  use var_inp,    only : infile, outfile
  use var_setp,   only : insimu
  use var_struct, only : nunit_all
  use var_mgo,    only : inmgo, iunit2sysmbr_mgo, iact2unit_mgo, &
                         ishadow2real_unit_mgo

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ----------------------------------------------------------------------
  ! intent(out) :: inmgo

  ! ----------------------------------------------------------------------
  ! function
  integer :: ifunc_char2int

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: j
  integer :: iunit, junit, icol, icol2
  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  integer :: isim, isys, istat, istat2, iactnum
  integer :: iact, jact, iact2sysmbr(3, MXACT_MGO)

  character(50) :: char50
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message


  ! ----------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  ! ----------------------------------------------------------------------
  ! default parameters
  inmgo%bdemax_mgo = -1.0
  inmgo%dihemax_mgo = -1.0
  inmgo%baemax_mgo = -1.0
  inmgo%enegap(1:MXSIM, 1:MXSYSTEM_MGO, 1:MXSTATE_MGO) = INVALID_VALUE
  inmgo%delta_mgo(1:MXSYSTEM_MGO, 1:MXSTATE_MGO, 1:MXSTATE_MGO) = INVALID_VALUE

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'multiple_go     ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "multiple_go" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     iunit = 0
     do iline = 1, nlines
        ctmp00 = cwkinp(iline)
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
   
        do iequa = 1, nequat

           cvalue = 'bdemax_mgo'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmgo%bdemax_mgo, cvalue)
   
           cvalue = 'baemax_mgo'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmgo%baemax_mgo, cvalue)
           
           cvalue = 'dihemax_mgo'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inmgo%dihemax_mgo, cvalue)

        end do
        
        if(ctmp00(1:6) == 'ENEGAP') then
           do icol = 7, CARRAY_MXCOLM
              if(ctmp00(icol:icol) == ')') exit
           end do
           read (ctmp00(8:icol-1), *) isim
           
           do icol2 = icol+2, CARRAY_MXCOLM
              if(ctmp00(icol2:icol2) == ')') exit
           end do
           read (ctmp00(icol+2:icol2-1), *) isys
           
           read (ctmp00, *) char50, &
                (inmgo%enegap(isim, isys, j), j = 1, inmgo%nstate_mgo(isys))
           write (lunout, '(2a)', advance='NO') '---reading ENEGAP: ', trim(char50)
           do j = 1, inmgo%nstate_mgo(isys)
              write (lunout, '(g10.3)', advance='NO') inmgo%enegap(isim, isys, j)
           end do
           write (lunout, *)

        else if(ctmp00(1:5) == 'DELTA') then
           do icol = 6, CARRAY_MXCOLM
              if(ctmp00(icol:icol) == ')') exit
           end do
           read (ctmp00(7:icol-3), *) isys
           istat = ifunc_char2int(ctmp00(icol-2:icol-2))
           istat2 = ifunc_char2int(ctmp00(icol-1:icol-1))

           read (ctmp00, *) char50, inmgo%delta_mgo(isys, istat, istat2)
           write (lunout, '(2a, g10.3)') '---reading DELTA: ', trim(char50), &
                inmgo%delta_mgo(isys, istat, istat2)

           if(inmgo%delta_mgo(isys, istat, istat2) >= 0.0) then
              inmgo%delta_mgo(isys, istat, istat2) = &
                   - inmgo%delta_mgo(isys, istat, istat2)
           end if
           inmgo%delta_mgo(isys, istat2, istat) = &
                inmgo%delta_mgo(isys, istat, istat2)
        end if
     end do


     ! -----------------------------------------------------------------
     ! checking input variables
     if(inmgo%bdemax_mgo <= 0.0) then
        error_message = 'Error: invalid value for bdemax_mgo'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inmgo%baemax_mgo <= 0.0) then
        error_message = 'Error: invalid value for baemax_mgo'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inmgo%dihemax_mgo <= 0.0) then
        error_message = 'Error: invalid value for dihemax_mgo'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do isim = 1, insimu%n_step_sim
        do isys = 1, inmgo%nsystem_mgo
           do istat = 1, inmgo%nstate_mgo(isys)
              if(inmgo%enegap(isim, isys, istat) > INVALID_JUDGE) then
                 error_message = 'Error: invalid value for enegap_mgo'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end do
        end do
     end do
     
     do isys = 1, inmgo%nsystem_mgo
        do istat = 1, inmgo%nstate_mgo(isys)
           do istat2 = istat + 1, inmgo%nstate_mgo(isys)
              if(inmgo%delta_mgo(isys, istat, istat2) > INVALID_JUDGE) then
                 error_message = 'Error: invalid value for delta_mgo'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end do
        end do
     end do



#ifdef MPI_PAR
  end if

!
! transfer size may be able to reduce. 
!
  call MPI_Bcast (inmgo, inmgo%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
#endif

  ! ----------------------------------------------------------------------
  ! calc iunit2sysmbr_mgo
  do iact = 1, MXACT_MGO
     iact2sysmbr(1, iact) = 0
     iact2sysmbr(2, iact) = 0
     iact2sysmbr(3, iact) = 0
  end do
  do isys = 1, inmgo%nsystem_mgo
     do istat = 1, inmgo%nstate_mgo(isys)
        do iactnum = 1, inmgo%nactnum_mgo(isys)
           iact = inmgo%isysmbr_mgo(isys, istat, iactnum)
           iact2sysmbr(1, iact) = isys
           iact2sysmbr(2, iact) = istat
           iact2sysmbr(3, iact) = iactnum
        end do
     end do
  end do
  
  do iunit = 1, nunit_all
     do junit = iunit, nunit_all
        iact = inmgo%iactmat_mgo(iunit, junit)
        if(iact /= 0) then
           isys = iact2sysmbr(1, iact)
           istat = iact2sysmbr(2, iact)
           iactnum = iact2sysmbr(3, iact)
           iunit2sysmbr_mgo(1, iunit, junit) = isys
           iunit2sysmbr_mgo(2, iunit, junit) = istat
           iunit2sysmbr_mgo(3, iunit, junit) = iactnum

        else
           iunit2sysmbr_mgo(1, iunit, junit) = 0 
           iunit2sysmbr_mgo(2, iunit, junit) = 0
           iunit2sysmbr_mgo(3, iunit, junit) = 0
       end if
     end do
  end do

  ! ----------------------------------------------------------------------
  ! making array: iact2unit_mgo
  do iunit = 1, nunit_all
     do junit = 1, nunit_all
        iact = inmgo%iactmat_mgo(iunit, junit)
        if(iact /= 0) then
           iact2unit_mgo(1, iact) = iunit
           iact2unit_mgo(2, iact) = junit
        end if
     end do
  end do

  ! ----------------------------------------------------------------------  
  ! making array: ishadow2real_unit_mgo
  do iunit = 1, nunit_all
     ishadow2real_unit_mgo(iunit) = iunit
  end do

  do isys = 1, inmgo%nsystem_mgo
     do istat = 2, inmgo%nstate_mgo(isys)
        do iactnum = 1, inmgo%nactnum_mgo(isys)
           iact = inmgo%isysmbr_mgo(isys, istat, iactnum)
           iunit = iact2unit_mgo(1, iact)
           junit = iact2unit_mgo(2, iact)

           if(iunit == junit) then
              jact = inmgo%isysmbr_mgo(isys, 1, iactnum)
              ishadow2real_unit_mgo(iunit) = iact2unit_mgo(1, jact)
           end if

        end do
     end do
  end do

end subroutine setp_para_mgo
