! mloop_setup_nlocal_mgo
!> @brief This subroutine is to perform the special treatment for the non-local interaction of multiple basins model.

! ***********************************************************************
subroutine mloop_setup_nlocal_mgo()

  use const_maxsize
  use const_index
  use var_inp,    only : outfile
  use var_struct, only : nmp_all, nunit_all, &
                         ncon, icon2unit, coef_go, ncon_unit
  use var_mgo,    only : inmgo, enegap_mgo, offset_mgo, offset_unit, &
                         !ireal2shadow_con_mgo, &
                         ncontype_mgo, icon2sysmbr_mgo, iunit2sysmbr_mgo
#ifdef MPI_PAR
  use mpiconst
#endif
  implicit none

  ! ----------------------------------------------------------------------
  ! intent(inout) :: ncon, icon2unit, ncon_unit, &
  ! go_nat, factor_go, offset_mgo, ncontype_mgo
  ! intent(out) :: ireal2shadow_con_mgo, icon2sysmbr_mgo

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: lunout
  integer :: iunit, junit, kcon, lcon
  integer :: isys, istat, iactnum, iact, jact, icon, l_ncon
  integer :: nact, jcon, iconnum
!  integer :: nact2con_mgo(MXACT_MGO)
!  integer :: lact2con_mgo(MXLCON_MGO, MXACT_MGO)
  integer :: lact2con_mgo(2, MXACT_MGO + 1)
  integer :: iact2con_mgo(nmp_all*MXMPCON)
  integer :: ipost2pre_con(nmp_all*MXMPCON), icopy(nmp_all*MXMPCON)
  character(CARRAY_MSG_ERROR) :: error_message
  
  ! ----------------------------------------------------------------------
  lunout = outfile%data

  ! ----------------------------------------------------------------------
  ! calc lact2con_mgo

  nact = inmgo%nact_mgo
  if (nact > MXACT_MGO) then
     error_message = 'Error: inmgo%nact_mgo > MXACT_MGO in mloop_setup_nlocal_mgo.'
     call util_error(ERROR%STOP_ALL, error_message)
  endif

  lact2con_mgo(1:2, 1:MXACT_MGO) = 0
  iconnum = 0
  do iactnum = 1, nact
     lact2con_mgo(1, iactnum) = iconnum + 1
     do icon = 1, ncon
        iunit = icon2unit(1, icon)
        junit = icon2unit(2, icon)
        iact = inmgo%iactmat_mgo(iunit, junit)
        if(iact == iactnum) then
           iconnum = iconnum + 1
           if(iconnum > nmp_all*MXMPCON) then
              error_message = 'Error: iconnum > nmp_all*MXMPCON in mloop_setup_nlocal_mgo.'
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           iact2con_mgo(iconnum) = icon
        end if
     end do
     lact2con_mgo(2, iactnum) = iconnum
  end do
  lact2con_mgo(1, nact + 1) = lact2con_mgo(2, nact) + 1
  lact2con_mgo(2, nact + 1) = lact2con_mgo(2, nact)

!  nact2con_mgo(1:MXACT_MGO) = 0
!  do icon = 1, ncon 
!     iunit = icon2unit(1, icon)
!     junit = icon2unit(2, icon)
!     iact = inmgo%iactmat_mgo(iunit, junit)
!     if(iact == 0) cycle
!     nact2con_mgo(iact) = nact2con_mgo(iact) + 1
!     if(nact2con_mgo(iact) > MXLCON_MGO) then
!        write (error_message, *) 'Error: too small MXLCON_MGO in mloop_setup_nlocal_mgo.'
!        call util_error(ERROR%STOP_ALL, error_message)
!     end if
!     lact2con_mgo(nact2con_mgo(iact), iact) = icon
!  end do

  ! ----------------------------------------------------------------------
  ! add contact
  l_ncon = ncon
  call mloop_setup_nlocal_add_mgo(l_ncon, lact2con_mgo, iact2con_mgo)
  ncon = l_ncon

  ! ----------------------------------------------------------------------
  ! sort contact
  call util_sort_contact(ipost2pre_con)
  do icon = 1, ncon
     icopy(icon) = ncontype_mgo(icon)
  end do

  do icon = 1, ncon
     ncontype_mgo(icon) = icopy(ipost2pre_con(icon))
  end do

  ! ----------------------------------------------------------------------
  ! recalc lact2con_mgo
  iconnum = 0
  do iactnum = 1, nact
     lact2con_mgo(1, iactnum) = iconnum + 1
     do icon = 1, ncon
        iunit = icon2unit(1, icon)
        junit = icon2unit(2, icon)
        iact = inmgo%iactmat_mgo(iunit, junit)
        if(iact == iactnum) then
           iconnum = iconnum + 1
           iact2con_mgo(iconnum) = icon
        end if
     end do
     lact2con_mgo(2, iactnum) = iconnum
  end do

!  nact2con_mgo(1:MXACT_MGO) = 0
!  do icon = 1, ncon 
!     iunit = icon2unit(1, icon)
!     junit = icon2unit(2, icon)
!     iact = inmgo%iactmat_mgo(iunit, junit)
!     if(iact == 0) cycle
!     nact2con_mgo(iact) = nact2con_mgo(iact) + 1
!     if(nact2con_mgo(iact) > MXLCON_MGO) then
!        error_message = 'Error: too small MXLCON_MGO in mloop_setup_nlocal_mgo'
!        call util_error(ERROR%STOP_ALL, error_message)
!     end if
!     lact2con_mgo(nact2con_mgo(iact), iact) = icon
!  end do

  ! ----------------------------------------------------------------------
  ! calc ncontype_mgo, irefcon_mgo
  call mloop_setup_nlocal_ref_mgo(lact2con_mgo, iact2con_mgo)

  ! ----------------------------------------------------------------------
  ! calc ireal2shadow_con_mgo
  !  !! If you want to use 'ireal2shadow_con_mgo', 
  !  !! comment-out next one line, and uncomment following do-enddo block.
  icon2sysmbr_mgo(1:2, 1:ncon) = 0
  !do icon = 1, nmp_all*MXMPCON
  !   icon2sysmbr_mgo(1, icon) = 0
  !   icon2sysmbr_mgo(2, icon) = 0
  !   ireal2shadow_con_mgo(1, icon) = 1
  !   do istat = 2, MXSTATE_MGO
  !      ireal2shadow_con_mgo(istat, icon) = 0
  !   end do
  !end do

  do isys = 1, inmgo%nsystem_mgo
     do iactnum = 1, inmgo%nactnum_mgo(isys)
        iact = inmgo%isysmbr_mgo(isys, 1, iactnum)

        do icon = lact2con_mgo(1, iact), lact2con_mgo(2, iact)
           kcon = iact2con_mgo(icon)
           icon2sysmbr_mgo(1, kcon) = isys
           icon2sysmbr_mgo(2, kcon) = 1
           !ireal2shadow_con_mgo(1, kcon) = inmgo%nstate_mgo(isys)     

           do istat = 2, inmgo%nstate_mgo(isys)                  
              jact = inmgo%isysmbr_mgo(isys, istat, iactnum)
              jcon = icon - lact2con_mgo(1, iact) + lact2con_mgo(1, jact)
              lcon = iact2con_mgo(jcon)
              icon2sysmbr_mgo(1, lcon) = isys
              icon2sysmbr_mgo(2, lcon) = istat
              !ireal2shadow_con_mgo(istat, kcon) = lcon
           end do

        end do
     end do
  end do

  ! ----------------------------------------------------------------------
  ! calc ncon_unit
  do iunit = 1, nunit_all
     do junit = 1, nunit_all
        ncon_unit(iunit, junit) = 0
     end do
  end do
  do icon = 1, ncon
     if(ncontype_mgo(icon) /= 2) then
        iunit = icon2unit(1, icon)
        junit = icon2unit(2, icon)
        ncon_unit(iunit, junit) = ncon_unit(iunit, junit) + 1 
     end if
  end do

  ! ----------------------------------------------------------------------
  ! For balancing the energy gap between different stats. 
  ! But following routine is appropriate,
  ! only when factor_go(icon) have same values. 
  ! ----------------------------------------------------------------------

  offset_mgo(1:inmgo%nsystem_mgo, 1:inmgo%nstate_max_mgo) = 0.0e0_PREC
  offset_unit(1:nunit_all, 1:nunit_all) = 0.0e0_PREC

  do icon = 1, ncon
     if (ncontype_mgo(icon) /= 2) then
        iunit = icon2unit(1, icon)
        junit = icon2unit(2, icon)
        iact  = inmgo%iactmat_mgo(iunit, junit)
        if (iact == 0) cycle
        isys  = iunit2sysmbr_mgo(1, iunit, junit)
        istat = iunit2sysmbr_mgo(2, iunit, junit)
        offset_mgo(isys, istat) = offset_mgo(isys, istat) + coef_go(icon)
        offset_unit(iunit, junit) = offset_unit(iunit, junit) + coef_go(icon)
     endif
  enddo

#ifdef MPI_PAR
  if (myrank == 0) then
#endif
  write (lunout, '(a)') '** for multiple Go,'
  write (lunout, '(a)') &
       '** energy gap between different states are offset. **'
  do isys = 1, inmgo%nsystem_mgo
     do istat = 1, inmgo%nstate_mgo(isys)
        write (lunout, '(a7, i2, a1, i2, a3, f10.3)') &
             'offset(', isys, ',', istat, ') =', offset_mgo(isys, istat), &
             'enegap(', isys, ',', istat, ') =', enegap_mgo(isys, istat)
     end do
  end do
#ifdef MPI_PAR
  endif
#endif

end subroutine mloop_setup_nlocal_mgo
