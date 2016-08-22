! mloop_setup_nlocal_add_mgo
!> @brief Create "dummy" contact in one state if the corresponding &
!>        contact naturally exists in other states. It is called by &
!>        mloop_setup_nlocal_mgo

subroutine mloop_setup_nlocal_add_mgo(l_ncon, lact2con_mgo, iact2con_mgo)

  use const_maxsize
  use const_index
  use var_struct, only : nmp_all, lunit2mp, icon2mp, icon2unit, ncon, &
                         factor_go, coef_go, icon_dummy_mgo
  use var_mgo,    only : inmgo, iact2unit_mgo, ncontype_mgo
  use var_setp,   only : inpro, inmisc
  implicit none

  ! ----------------------------------------------------------------------
  integer, intent(inout) :: l_ncon, lact2con_mgo(2, MXACT_MGO + 1)
  integer, intent(inout) :: iact2con_mgo(nmp_all*MXMPCON)
  ! intent(inout) :: icon2mp, icon2unit, factor_go
  ! intent(inout) :: ncontype(MXCON)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: ii, jj
  integer :: imp1, imp2, jmp1, jmp2
  integer :: isys, istat, jstat, iactnum
  integer :: icon, jcon, iunit, junit, nact
  integer :: iact(MXSTATE_MGO), addcon(MXACT_MGO)
  integer :: iunit1(MXSTATE_MGO), iunit2(MXSTATE_MGO)
  integer :: ipre1(MXSTATE_MGO), ipre2(MXSTATE_MGO)
  character(CARRAY_MSG_ERROR) :: error_message

  ! ----------------------------------------------------------------------
  nact = inmgo%nact_mgo
  addcon(1:MXACT_MGO) = 0
  do icon = 1, ncon
     if (icon_dummy_mgo(icon) == 0) then
        ncontype_mgo(icon) = 2 
     else 
        ncontype_mgo(icon) = 0
     endif
  enddo

  ! ----------------------------------------------------------------------
  ! add contact due to the multiple Go potential
  do isys = 1, inmgo%nsystem_mgo
     if(inmgo%nstate_mgo(isys) == 1) cycle
     do iactnum = 1, inmgo%nactnum_mgo(isys)

        ! ----------------------------------------------------------------
        do istat = 1, inmgo%nstate_mgo(isys)
           iact(istat) = inmgo%isysmbr_mgo(isys, istat, iactnum)
           iunit1(istat) = iact2unit_mgo(1, iact(istat)) 
           iunit2(istat) = iact2unit_mgo(2, iact(istat))
           ipre1(istat) = lunit2mp(1, iunit1(istat)) - 1
           ipre2(istat) = lunit2mp(1, iunit2(istat)) - 1
        end do

        do istat = 1, inmgo%nstate_mgo(isys)
           do ii = lact2con_mgo(1, iact(istat)), lact2con_mgo(2, iact(istat))
              icon = iact2con_mgo(ii)
              imp1 = icon2mp(1, icon) - ipre1(istat)
              imp2 = icon2mp(2, icon) - ipre2(istat)
              loop_jstat: do jstat = 1, inmgo%nstate_mgo(isys)
                 if(jstat == istat) cycle loop_jstat
                 do jj = lact2con_mgo(1, iact(jstat)), lact2con_mgo(2, iact(jstat))
                    jcon = iact2con_mgo(jj) 
                    jmp1 = icon2mp(1, jcon) - ipre1(jstat)
                    jmp2 = icon2mp(2, jcon) - ipre2(jstat)
                    if(imp1 == jmp1 .and. imp2 == jmp2) cycle loop_jstat
                 end do
                 do jj = lact2con_mgo(1, nact + 1), lact2con_mgo(2, nact + 1)
                    jcon = iact2con_mgo(jj) 
                    iunit = icon2unit(1, jcon)
                    junit = icon2unit(2, jcon)
                    if(iunit == iunit1(jstat) .and. junit == iunit2(jstat)) then
                       jmp1 = icon2mp(1, jcon) - ipre1(jstat)
                       jmp2 = icon2mp(2, jcon) - ipre2(jstat)
                       if(imp1 == jmp1 .and. imp2 == jmp2) cycle loop_jstat
                    end if
                 end do

                 ! add contact
                 l_ncon = l_ncon + 1
                 icon2mp(1, l_ncon) = imp1 + ipre1(jstat)
                 icon2mp(2, l_ncon) = imp2 + ipre2(jstat)
                 icon2unit(1, l_ncon) = iunit1(jstat)
                 icon2unit(2, l_ncon) = iunit2(jstat)
                 icon_dummy_mgo(l_ncon) = 0
                 factor_go(l_ncon) = inmisc%factor_go_unit(iunit1(jstat), iunit2(jstat))
                 coef_go(l_ncon) = factor_go(l_ncon) * inpro%cgo1210
                 ncontype_mgo(l_ncon) = 2
                 lact2con_mgo(2, nact + 1) = lact2con_mgo(2, nact + 1) + 1
                 iact2con_mgo(lact2con_mgo(2, nact + 1)) = l_ncon
                 addcon(iact(jstat)) = addcon(iact(jstat)) + 1
              end do loop_jstat
           end do
        end do

        ! ----------------------------------------------------------------
        ! check the number of contact for each act
        do istat = 1, inmgo%nstate_mgo(isys)
           do jstat = istat + 1, inmgo%nstate_mgo(isys)
              if(lact2con_mgo(2, iact(istat)) - lact2con_mgo(1, iact(istat)) + addcon(iact(istat)) &
                   /= lact2con_mgo(2, iact(jstat)) - lact2con_mgo(1, iact(jstat)) + addcon(iact(jstat))) then
                 error_message = 'Error: different contact number between state in mloop_setup_nlocal_add_mgo'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end do
        end do

     end do
  end do

end subroutine mloop_setup_nlocal_add_mgo
