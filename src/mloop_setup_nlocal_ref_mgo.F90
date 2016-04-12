! mloop_setup_nlocal_ref_mgo
!> @brief Make Go potential for multiple Go model

! ***********************************************************************
subroutine mloop_setup_nlocal_ref_mgo(lact2con_mgo, iact2con_mgo)

  use const_maxsize
  use var_struct, only : ncon, lunit2mp, icon2mp, icon2unit, go_nat, go_nat2
  use var_mgo, only : inmgo, iact2unit_mgo, ncontype_mgo, irefcon_mgo
  implicit none

  ! ----------------------------------------------------------------------
  integer, intent(in) :: lact2con_mgo(2, MXACT_MGO + 1)
  integer, intent(in) :: iact2con_mgo(MXCON)
  ! intent(inout) :: icon2mp, icon2unit, ncontype_mgo
  ! intent(out) :: irefcon_mgo

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: ii
  integer :: isys, istat, iactnum, kcon
  integer :: iact(MXSTATE_MGO), icon(MXSTATE_MGO)
  integer :: iunit1(MXSTATE_MGO), iunit2(MXSTATE_MGO)
  integer :: ipre1(MXSTATE_MGO), ipre2(MXSTATE_MGO)
  integer :: imingo_stat
  real(PREC) :: mingo_nat

  ! ----------------------------------------------------------------------
  ! calc ncontype_mgo and irefcon_mgo
  do kcon = 1, ncon
     irefcon_mgo(kcon) = 0
  end do

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

        do ii = lact2con_mgo(1, iact(1)), lact2con_mgo(2, iact(1))
           mingo_nat = 1000.0e0_PREC
           imingo_stat = 0
           do istat = 1, inmgo%nstate_mgo(isys)
              icon(istat) = iact2con_mgo(ii- lact2con_mgo(1, iact(1)) + lact2con_mgo(1, iact(istat)))
              if(ncontype_mgo(icon(istat)) == 2) cycle
              if(go_nat(icon(istat)) < mingo_nat) then
                 mingo_nat = go_nat(icon(istat))
                 imingo_stat = istat
              end if
           end do

           do istat = 1, inmgo%nstate_mgo(isys)
              ! not contact
              if(ncontype_mgo(icon(istat)) == 2) then
                 ncontype_mgo(icon(istat)) = 2
                 go_nat(icon(istat)) = mingo_nat
                 go_nat2(icon(istat)) = go_nat(icon(istat))**2
              ! the nearest contact
              else if(istat == imingo_stat) then
                 ncontype_mgo(icon(istat)) = 0
              ! not the nearest contact
              else
                 ncontype_mgo(icon(istat)) = 1
              end if
              irefcon_mgo(icon(istat)) = icon(imingo_stat)
           end do
        end do

     end do
  end do

end subroutine mloop_setup_nlocal_ref_mgo
