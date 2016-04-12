! mloop_setup_local_mgo
!> @brief Setting local interactions for multiple Go model

! ***********************************************************************         
subroutine mloop_setup_local_mgo()

  use const_maxsize
  use const_index
  use var_struct, only : nmp_all, imp2unit, iclass_mp, &
                         nbd,  ibd2mp,  bd_nat,  coef_bd,  correct_bd_mgo, &
                         nba,  iba2mp,  ba_nat,  coef_ba,  correct_ba_mgo, &
                         ndih, idih2mp, dih_nat, coef_dih, correct_dih_mgo
  use var_mgo, only : inmgo,          & 
                      ibd2sysmbr_mgo, &
                      iba2sysmbr_mgo, &
                      idih2sysmbr_mgo
  implicit none

  ! ----------------------------------------------------------------------
  ! intent(inout) :: factor_bd, factor_ba, factor_dih

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: kmp1, kmp2, inum, jnum
  integer :: kbd, kba, kdih, kact, kunit1, kunit2
  integer :: isys, istat, jstat, iactnum
  integer :: iact(MXSTATE_MGO)
  integer :: nact, ibdnum, ibanum, idihnum
!  integer :: nact2bd_mgo(MXACT_MGO), lact2bd_mgo(MXBD, MXACT_MGO)
!  integer :: nact2ba_mgo(MXACT_MGO), lact2ba_mgo(MXBA, MXACT_MGO)
!  integer :: nact2dih_mgo(MXACT_MGO), lact2dih_mgo(MXDIH, MXACT_MGO)
  integer :: lact2bd_mgo(2, MXACT_MGO),  iact2bd_mgo(MXMPBD*nmp_all)
  integer :: lact2ba_mgo(2, MXACT_MGO),  iact2ba_mgo(MXMPBA*nmp_all)
  integer :: lact2dih_mgo(2, MXACT_MGO), iact2dih_mgo(MXMPDIH*nmp_all)
  real(PREC) :: e1, maxe1, eratio, ddist2

  ! ----------------------------------------------------------------------
  correct_bd_mgo(1:MXMPBD*nmp_all) = 1.0
  correct_ba_mgo(1:MXMPBA*nmp_all) = 1.0
  correct_dih_mgo(1:MXMPDIH*nmp_all) = 1.0

  nact = inmgo%nact_mgo

  lact2bd_mgo(1:2, 1:MXACT_MGO) = 0
  ibdnum = 0
  lact2ba_mgo(1:2, 1:MXACT_MGO) = 0
  ibanum = 0
  lact2dih_mgo(1:2, 1:MXACT_MGO) = 0
  idihnum = 0
  ibd2sysmbr_mgo(1:2, 1:nbd) = 0
  iba2sysmbr_mgo(1:2, 1:nba) = 0
  idih2sysmbr_mgo(1:2, 1:ndih) = 0
  do iactnum = 1, nact
     !--  bond  -----------------------------
     lact2bd_mgo(1, iactnum) = ibdnum + 1
     do kbd = 1, nbd
        kmp1 = ibd2mp(1, kbd)
        kmp2 = ibd2mp(2, kbd)
        kunit1 = imp2unit(kmp1)
        kunit2 = imp2unit(kmp2)
        kact = inmgo%iactmat_mgo(kunit1, kunit2)
        if (kact == iactnum) then
           ibdnum = ibdnum + 1
           iact2bd_mgo(ibdnum) = kbd
        endif
     end do
     lact2bd_mgo(2, iactnum) = ibdnum

     !--  bond-angle  -----------------------------
     lact2ba_mgo(1, iactnum) = ibanum + 1
     do kba = 1, nba
        kmp1 = iba2mp(1, kba)
        kmp2 = iba2mp(3, kba)
        kunit1 = imp2unit(kmp1)
        kunit2 = imp2unit(kmp2)
        kact = inmgo%iactmat_mgo(kunit1, kunit2)
        if (kact == iactnum) then
           ibanum = ibanum + 1
           iact2ba_mgo(ibanum) = kba
        endif
     end do
     lact2ba_mgo(2, iactnum) = ibanum

     !--  dihedral-angle  -----------------------------
     lact2dih_mgo(1, iactnum) = idihnum + 1
     do kdih = 1, ndih
        kmp1 = idih2mp(1, kdih)
        kmp2 = idih2mp(4, kdih)
        kunit1 = imp2unit(kmp1)
        kunit2 = imp2unit(kmp2)
        kact = inmgo%iactmat_mgo(kunit1, kunit2)
        if (kact == iactnum) then
           idihnum = idihnum + 1
           iact2dih_mgo(idihnum) = kdih
        endif
     end do
     lact2dih_mgo(2, iactnum) = idihnum

  enddo

!  do kact = 1, MXACT_MGO
!     nact2bd_mgo(kact) = 0
!  end do
!  do kbd = 1, nbd
!     ibd2sysmbr_mgo(1, kbd) = 0
!     ibd2sysmbr_mgo(2, kbd) = 0
!     kmp1 = ibd2mp(1, kbd)
!     kmp2 = ibd2mp(2, kbd)
!     kunit1 = imp2unit(kmp1)
!     kunit2 = imp2unit(kmp2)
!     kact = inmgo%iactmat_mgo(kunit1, kunit2)
!     if(kact == 0) cycle
!     nact2bd_mgo(kact) = nact2bd_mgo(kact) + 1
!     lact2bd_mgo(nact2bd_mgo(kact), kact) = kbd
!  end do
!
!  do kact = 1, MXACT_MGO
!     nact2ba_mgo(kact) = 0
!  end do
!  do kba = 1, nba
!     iba2sysmbr_mgo(1, kba) = 0
!     iba2sysmbr_mgo(2, kba) = 0
!     kmp1 = iba2mp(1, kba)
!     kmp2 = iba2mp(3, kba)
!     kunit1 = imp2unit(kmp1)
!     kunit2 = imp2unit(kmp2)
!     kact = inmgo%iactmat_mgo(kunit1, kunit2)
!     if(kact == 0) cycle
!     nact2ba_mgo(kact) = nact2ba_mgo(kact) + 1
!     lact2ba_mgo(nact2ba_mgo(kact), kact) = kba
!  end do
!
!  do kact = 1, MXACT_MGO
!     nact2dih_mgo(kact) = 0
!  end do
!  do kdih = 1, ndih
!     idih2sysmbr_mgo(1, kdih) = 0
!     idih2sysmbr_mgo(2, kdih) = 0
!     kmp1 = idih2mp(1, kdih)
!     kmp2 = idih2mp(4, kdih)
!     kunit1 = imp2unit(kmp1)
!     kunit2 = imp2unit(kmp2)
!     kact = inmgo%iactmat_mgo(kunit1, kunit2)
!     if(kact == 0) cycle
!     nact2dih_mgo(kact) = nact2dih_mgo(kact) + 1
!     lact2dih_mgo(nact2dih_mgo(kact), kact) = kdih
!  end do

  ! ----------------------------------------------------------------------
  ! modify factor_bd, factor_ba, factor_dih
  do isys = 1, inmgo%nsystem_mgo
     do iactnum = 1, inmgo%nactnum_mgo(isys)
        do istat = 1, inmgo%nstate_mgo(isys)
           iact(istat) = inmgo%isysmbr_mgo(isys, istat, iactnum)
        end do

        !--  bond  -----------------------------
        !do kbd = 1, nact2bd_mgo(iact(1))
        do kbd = lact2bd_mgo(1, iact(1)), lact2bd_mgo(2, iact(1))
           maxe1 = 0.0e0_PREC
           do istat = 1, inmgo%nstate_mgo(isys)
              !inum = lact2bd_mgo(kbd, iact(istat))
              inum = iact2bd_mgo(kbd-lact2bd_mgo(1,iact(1))+lact2bd_mgo(1,iact(istat)))
              do jstat = istat + 1, inmgo%nstate_mgo(isys)
                 !jnum = lact2bd_mgo(kbd, iact(jstat))
                 jnum = iact2bd_mgo(kbd-lact2bd_mgo(1,iact(1))+lact2bd_mgo(1,iact(jstat)))
                 ddist2 = (bd_nat(jnum) - bd_nat(inum))**2
                 e1 = (coef_bd(1, inum) + coef_bd(2, inum) * ddist2) * ddist2
                 if(e1 > maxe1) maxe1 = e1
              end do
              ibd2sysmbr_mgo(1, inum) = isys
              ibd2sysmbr_mgo(2, inum) = istat
           end do
           if(maxe1 > inmgo%bdemax_mgo) then
              eratio = inmgo%bdemax_mgo / maxe1
              do istat = 1, inmgo%nstate_mgo(isys)
                 !inum = lact2bd_mgo(kbd, iact(istat))
                 inum = iact2bd_mgo(kbd-lact2bd_mgo(1,iact(1))+lact2bd_mgo(1,iact(istat)))
                 correct_bd_mgo(inum) = eratio
                 coef_bd(1, inum) = coef_bd(1, inum) * eratio
                 coef_bd(2, inum) = coef_bd(2, inum) * eratio
              end do
           end if
        end do

        !--  bond-angle  -----------------------------
        !do kba = 1, nact2ba_mgo(iact(1))
        do kba = lact2ba_mgo(1, iact(1)), lact2ba_mgo(2, iact(1))
           maxe1 = 0.0e0_PREC
           do istat = 1, inmgo%nstate_mgo(isys)
              !inum = lact2ba_mgo(kba, iact(istat))
              inum = iact2ba_mgo(kba-lact2ba_mgo(1,iact(1))+lact2ba_mgo(1,iact(istat)))
              do jstat = istat + 1, inmgo%nstate_mgo(isys)
                 !jnum = lact2ba_mgo(kba, iact(jstat))
                 jnum = iact2ba_mgo(kba-lact2ba_mgo(1,iact(1))+lact2ba_mgo(1,iact(jstat)))
                 ddist2 = (ba_nat(jnum) - ba_nat(inum))**2
                 e1 = coef_ba(1, inum) * ddist2
                 if(e1 > maxe1) maxe1 = e1
              end do
              iba2sysmbr_mgo(1, inum) = isys
              iba2sysmbr_mgo(2, inum) = istat
           end do
           if(maxe1 > inmgo%baemax_mgo) then
              eratio = inmgo%baemax_mgo / maxe1
              do istat = 1, inmgo%nstate_mgo(isys)
                 !inum = lact2ba_mgo(kba, iact(istat))
                 inum = iact2ba_mgo(kba-lact2ba_mgo(1,iact(1))+lact2ba_mgo(1,iact(istat)))
                 correct_ba_mgo(inum) = eratio
                 coef_ba(1, inum) = coef_ba(1, inum) * eratio
                 coef_ba(2, inum) = coef_ba(2, inum) * eratio
              end do
           end if
        end do

        !--  dihedral-angle  -----------------------------
        !do kdih = 1, nact2dih_mgo(iact(1))
        do kdih = lact2dih_mgo(1, iact(1)), lact2dih_mgo(2,iact(1))
           inum = iact2dih_mgo(kdih)
           kmp1 = idih2mp(1, inum)
           kmp2 = idih2mp(4, inum)
           if(iclass_mp(kmp1) == CLASS%LIG .AND. &
              iclass_mp(kmp2) == CLASS%LIG) cycle

           maxe1 = 0.0e0_PREC
           do istat = 1, inmgo%nstate_mgo(isys)
              !inum = lact2dih_mgo(kdih, iact(istat))
              inum = iact2dih_mgo(kdih-lact2dih_mgo(1,iact(1))+lact2dih_mgo(1,iact(istat)))
              do jstat = istat + 1, inmgo%nstate_mgo(isys)
                 !jnum = lact2dih_mgo(kdih, iact(jstat))
                 jnum = iact2dih_mgo(kdih-lact2dih_mgo(1,iact(1))+lact2dih_mgo(1,iact(jstat)))
                 ddist2 = (dih_nat(jnum) - dih_nat(inum))
                 e1 = coef_dih(1, inum) * (1 - cos(ddist2)) + &
                      coef_dih(2, inum) * (1 - cos(3*ddist2))
                 if(e1 > maxe1) maxe1 = e1
              end do
              idih2sysmbr_mgo(1, inum) = isys
              idih2sysmbr_mgo(2, inum) = istat
           end do
           if(maxe1 > inmgo%dihemax_mgo) then
              eratio = inmgo%dihemax_mgo / maxe1
              do istat = 1, inmgo%nstate_mgo(isys)
                 !inum = lact2dih_mgo(kdih, iact(istat))
                 inum = iact2dih_mgo(kdih-lact2dih_mgo(1,iact(1))+lact2dih_mgo(1,iact(istat)))
                 correct_dih_mgo(inum) = eratio
                 coef_dih(1, inum) = coef_dih(1, inum) * eratio
                 coef_dih(2, inum) = coef_dih(2, inum) * eratio
              end do
           end if
        end do

        !-- harmonic dihedral-angle  -----------------------------
        !do kdih = 1, nact2dih_mgo(iact(1))
        do kdih = lact2dih_mgo(1, iact(1)), lact2dih_mgo(2,iact(1))
           inum = iact2dih_mgo(kdih)
           kmp1 = idih2mp(1, inum)
           kmp2 = idih2mp(4, inum)
           if(iclass_mp(kmp1) /= CLASS%LIG .OR. &
              iclass_mp(kmp2) /= CLASS%LIG) cycle

           maxe1 = 0.0e0_PREC
           do istat = 1, inmgo%nstate_mgo(isys)
              !inum = lact2dih_mgo(kdih, iact(istat))
              inum = iact2dih_mgo(kdih-lact2dih_mgo(1,iact(1))+lact2dih_mgo(1,iact(istat)))
              do jstat = istat + 1, inmgo%nstate_mgo(isys)
                 !jnum = lact2dih_mgo(kdih, iact(jstat))
                 jnum = iact2dih_mgo(kdih-lact2dih_mgo(1,iact(1))+lact2dih_mgo(1,iact(jstat)))
                 ddist2 = (acos(cos(dih_nat(jnum))) - acos(cos(dih_nat(inum))))**2
                 e1 = coef_dih(1, inum) * ddist2
                 if(e1 > maxe1) maxe1 = e1
              end do
              idih2sysmbr_mgo(1, inum) = isys
              idih2sysmbr_mgo(2, inum) = istat
           end do
           if(maxe1 > inmgo%dihemax_mgo) then
              eratio = inmgo%dihemax_mgo / maxe1
              do istat = 1, inmgo%nstate_mgo(isys)
                 !inum = lact2dih_mgo(kdih, iact(istat))
                 inum = iact2dih_mgo(kdih-lact2dih_mgo(1,iact(1))+lact2dih_mgo(1,iact(istat)))
                 correct_dih_mgo(inum) = eratio
                 coef_dih(1, inum) = coef_dih(1, inum) * eratio
                 coef_dih(2, inum) = 0.0e0_PREC
              end do
           end if
        end do
        
     end do
  end do

end subroutine mloop_setup_local_mgo
