! read_nativeinfo
!> @brief Reading nativeinfo file

! **********************************************************************
! This subroutine is for reading intra coordinates from the
! file defined by lun
! **********************************************************************
subroutine read_nativeinfo(lun, i_ninfo_type, iunit, junit)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp, only : inmisc
  use var_struct, only : lunit2mp, imp2unit, cmp2seq,&
       nbd, ibd2mp, bd_nat, factor_bd, coef_bd, &
       nfene, ifene2mp, fene_nat, coef_fene, dist2_fene, &
       nrouse, irouse2mp, coef_rouse, &
       nba, iba2mp, ba_nat, factor_ba, coef_ba, &
!       ndih, idih2mp, dih_nat, dih_sin_nat, dih_cos_nat, factor_dih, coef_dih, idih2type, &
       ncon, icon2mp, icon2unit, go_nat, go_nat2, factor_go, icon_dummy_mgo, coef_go, &
       nLJ, iLJ2mp, iLJ2unit, LJ_nat, LJ_nat2, coef_LJ,&
       nwca, iwca2mp, iwca2unit, wca_nat, wca_nat2, coef_wca,&
       ncon_gauss, icon_gauss2mp, icon_gauss2unit, &
!       nrna_bp, irna_bp2mp, irna_bp2unit, rna_bp_nat, rna_bp_nat2, &
!       factor_rna_bp, irna_bp_dummy_mgo, nhb_bp, coef_rna_bp, &
!       nrna_st, irna_st2mp, irna_st2unit, rna_st_nat, rna_st_nat2, &
!       factor_rna_st, irna_st_dummy_mgo, coef_rna_st, &
       ibd2type, iba2type, icon2type, &
       idtrna_st2mp, dtrna_st_nat, coef_dtrna_st, ndtrna_st, idtrna_st2nn, &
       idtrna_hb2mp, dtrna_hb_nat, coef_dtrna_hb, ndtrna_hb, &
       idtrna_hb2hbsite, imp2hbsite, flg_hb_tertiary,&
       ndtrna_tst, idtrna_tst2mp, dtrna_tst_nat, coef_dtrna_tst,&
       flg_tst_exclusive, idtrna_tst2side

  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lun, i_ninfo_type, iunit, junit
  ! intent(out) :: nbd, ibd2mp, bd_nat, factor_bd, &
  !     nba, iba2mp, ba_nat, factor_ba, &
  !     ndih, idih2mp, dih_nat, factor_dih, &
  !     ncon, icon2mp, icon2unit, go_nat, factor_go

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: input_status
  integer :: ii, jj, imp1, imp2, imp3, imp4, iunit1, iunit2, imp_tmp
  integer :: imp1un, imp2un, imp3un, imp4un, kunit1
  integer :: ibd, ifene, irouse, iba, icon, iLJ, iwca, icon_gauss
!  integer :: idih
  integer :: ibd_read, iba_read, idih_read
  integer :: icon_read, idummy
!  integer :: ibp, ibp_read, ist
  integer :: nHB, ist_read
  integer :: ibsdist
  integer :: ihbdist, ihb_read
  integer :: itbsdist
!  integer :: itype
  integer :: offset1, offset2
  integer :: ist_type1, ist_type2, ist_type_tmp
  real(PREC) :: bl, ba, dih, go, factor, correct, coef !, coef3
  real(PREC) :: dist, energy0
  character(256) :: cline, cline_head
  character(CARRAY_MSG_ERROR) :: error_message
  character(1) :: ctype1
  character(2) :: ctype2
  character(3) :: ctype3
  character(4) :: ctype4
  character(4) :: a11,a12,a21,a22,a31,a32, atmp

  ! For checking if all angl and dihd are stored for hb-dist, st-dist, tst-dist
  logical :: datacheck_dtrna_hb_angl(2,MXDTRNAHB)
  logical :: datacheck_dtrna_hb_dihd(3,MXDTRNAHB)
  logical :: datacheck_dtrna_st_dihd(2,MXDTRNAST)
  logical :: datacheck_dtrna_tst_angl(2,MXDTRNATST)
  logical :: datacheck_dtrna_tst_dihd(3,MXDTRNATST)

  integer ::ifunc_nn2id
  ! ---------------------------------------------------------------------

  ibd = nbd 
  ifene = nfene
  irouse = nrouse
  iba = nba 
!  idih = ndih
  icon = ncon
  iLJ = nLJ
  iwca = nwca
  icon_gauss = ncon_gauss
!  ibp = nrna_bp
!  ist = nrna_st
  ibsdist = ndtrna_st
  ihbdist = ndtrna_hb
  itbsdist = ndtrna_tst
  ii = lunit2mp(1, iunit) - 1
  jj = lunit2mp(1, junit) - 1

  !! If XXX-diste is read, corresponding position of this flag is changed to False
  !! Then once angl/dihd for that dist is read, the flag is changed back to True.
  !! At last, all elements have to be True.
  datacheck_dtrna_hb_angl(1:2, 1:MXDTRNAHB) = .True.
  datacheck_dtrna_hb_dihd(1:3, 1:MXDTRNAHB) = .True.
  datacheck_dtrna_st_dihd(1:2, 1:MXDTRNAST) = .True.
  datacheck_dtrna_tst_angl(1:2, 1:MXDTRNATST) = .True.
  datacheck_dtrna_tst_dihd(1:3, 1:MXDTRNATST) = .True.

  ! ---------------------------------------------------------------------
  ! reading intra coodinates
  ! ---------------------------------------------------------------------
  do
     read (lun, '(a256)', iostat = input_status) cline
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: cannot read intra coordinate in read_nativeinfo'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     ! ------------------------------------------------------------------
     ! read the bond length
     if(cline(1:4) == 'bond') then
        ctype2 = '  '
        read (cline, *, iostat = input_status)     &
             cline_head, ibd_read, iunit1, iunit2, &
             imp1, imp2, imp1un, imp2un,           &
             bl, factor, correct, coef, ctype2
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        ibd = ibd + 1
        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           ibd2mp(1, ibd) = imp1 + ii
           ibd2mp(2, ibd) = imp2 + ii
           kunit1 = iunit1
        else ! NATIVEINFO%ONE_BY_ONE
           ibd2mp(1, ibd) = imp1un + ii
           ibd2mp(2, ibd) = imp2un + ii
           kunit1 = iunit
        end if
        bd_nat(ibd) = bl 
        factor_bd(ibd) = factor
        if (inmisc%flg_coef_from_ninfo) then
           coef_bd(1,ibd ) = coef
        endif

        if (ctype2 /= '  ') then
           ibd2type(ibd) = str2bondtype(ctype2)
        endif
     end if

     ! ------------------------------------------------------------------
     ! read FENE
     if(cline(1:4) == 'fene') then
        ctype2 = '  '
        read (cline, *, iostat = input_status)     &
             cline_head, ibd_read, iunit1, iunit2, &
             imp1, imp2, imp1un, imp2un,           &
             bl, factor, coef, ctype2
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        ifene = ifene + 1
        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           ifene2mp(1, ifene) = imp1 + ii
           ifene2mp(2, ifene) = imp2 + ii
           kunit1 = iunit1
        else ! NATIVEINFO%ONE_BY_ONE
           ifene2mp(1, ifene) = imp1un + ii
           ifene2mp(2, ifene) = imp2un + ii
           kunit1 = iunit
        end if
        fene_nat(ifene) = bl 
        if (inmisc%flg_coef_from_ninfo) then
           coef_fene(ifene) = coef
           dist2_fene(ifene) = factor
        endif

        !if (ctype2 /= '  ') then
        !   ibd2type(ibd) = str2bondtype(ctype2)
        !endif
     end if

     ! ------------------------------------------------------------------
     ! read the Rouse spring (bond)
     if(cline(1:5) == 'rouse') then
        read (cline, *, iostat = input_status)     &
             cline_head, ibd_read, iunit1, iunit2, &
             imp1, imp2, imp1un, imp2un, coef
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        irouse = irouse + 1
        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           irouse2mp(1, irouse) = imp1 + ii
           irouse2mp(2, irouse) = imp2 + ii
           kunit1 = iunit1
        else ! NATIVEINFO%ONE_BY_ONE
           irouse2mp(1, irouse) = imp1un + ii
           irouse2mp(2, irouse) = imp2un + ii
           kunit1 = iunit
        end if
        !if (inmisc%flg_coef_from_ninfo) then
           coef_rouse(1,irouse,:) = coef
        !endif

     end if

     ! ------------------------------------------------------------------
     ! read the bond angle
     if(cline(1:4) == 'angl') then
        ctype3 = '   '
        read (cline, *, iostat = input_status) &
            cline_head, iba_read, iunit1, iunit2,     &
            imp1, imp2, imp3, imp1un, imp2un, imp3un, &
            ba, factor, correct, coef, ctype3
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        iba = iba + 1 
        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           iba2mp(1, iba) = imp1 + ii
           iba2mp(2, iba) = imp2 + ii
           iba2mp(3, iba) = imp3 + ii
        else ! NATIVEINFO%ONE_BY_ONE
           iba2mp(1, iba) = imp1un + ii
           iba2mp(2, iba) = imp2un + ii
           iba2mp(3, iba) = imp3un + ii
        end if
        ba_nat(iba) = (ba / 180.0e0_PREC) * F_PI
        factor_ba(iba) = factor
        if (inmisc%flg_coef_from_ninfo) then
           coef_ba(1,iba) = coef
        endif

        if (ctype3 /= '   ') then
           iba2type(iba) = str2angletype(ctype3)
        endif
     end if

     ! ------------------------------------------------------------------
     ! read contact
     if(cline(1:7) == 'contact') then
        ctype3 = '   '
        read (cline, *, iostat = input_status) &
             cline_head, icon_read, iunit1, iunit2,  &
             imp1, imp2, imp1un, imp2un,             &
             go, factor, idummy, coef, ctype3
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        icon = icon + 1
        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + jj
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + jj
        end if
        icon2unit(1, icon) = imp2unit(imp1)
        icon2unit(2, icon) = imp2unit(imp2)
        icon2mp(1, icon) = imp1
        icon2mp(2, icon) = imp2
        go_nat(icon) = go
        go_nat2(icon) = go**2
        factor_go(icon) = factor
        icon_dummy_mgo(icon) = idummy
        if (inmisc%flg_coef_from_ninfo) then
           coef_go(icon) = coef
        endif

        if (ctype3 /= '   ') then
           icon2type(icon) = str2gotype(ctype3)
        endif
     end if

     ! ------------------------------------------------------------------
     ! read LJ
     if(cline(1:2) == 'LJ') then
        read (cline, *, iostat = input_status) &
             cline_head, icon_read, iunit1, iunit2,  &
             imp1, imp2, imp1un, imp2un,             &
             go, coef
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        iLJ = iLJ + 1
        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + jj
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + jj
        end if
        iLJ2unit(1, iLJ) = imp2unit(imp1)
        iLJ2unit(2, iLJ) = imp2unit(imp2)
        iLJ2mp(1, iLJ) = imp1
        iLJ2mp(2, iLJ) = imp2
        LJ_nat(iLJ) = go
        LJ_nat2(iLJ) = go**2
        if (inmisc%flg_coef_from_ninfo) then
           coef_LJ(iLJ) = coef
        endif
     end if

     ! ------------------------------------------------------------------
     ! read wca
     if(cline(1:3) == 'wca') then
        read (cline, *, iostat = input_status) &
             cline_head, icon_read, iunit1, iunit2,  &
             imp1, imp2, imp1un, imp2un,             &
             go, coef, factor
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        iwca = iwca + 1
        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + jj
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + jj
        end if
        iwca2unit(1, iwca) = imp2unit(imp1)
        iwca2unit(2, iwca) = imp2unit(imp2)
        iwca2mp(1, iwca) = imp1
        iwca2mp(2, iwca) = imp2
        wca_nat(iwca) = go
        wca_nat2(iwca) = go**2
        !if (inmisc%flg_coef_from_ninfo) then
           coef_wca(iwca,1) = coef
           coef_wca(iwca,2) = factor
        !endif
     end if

     ! ------------------------------------------------------------------
     ! read con_gauss
     if(cline(1:9) == 'con_gauss') then
        read (cline, *, iostat = input_status) &
             cline_head, icon_read, iunit1, iunit2,  &
             imp1, imp2, imp1un, imp2un
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        icon_gauss = icon_gauss + 1
        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + jj
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + jj
        end if
        icon_gauss2unit(1, icon_gauss) = imp2unit(imp1)
        icon_gauss2unit(2, icon_gauss) = imp2unit(imp2)
        icon_gauss2mp(1, icon_gauss) = imp1
        icon_gauss2mp(2, icon_gauss) = imp2
        !if (inmisc%flg_coef_from_ninfo) then
           !coef_con_gauss(icon_gauss) = coef
        !endif
     end if

!     ! ------------------------------------------------------------------
!     ! read basepair
!     if(cline(1:8) == 'basepair') then
!        ctype3 = '   '
!        nHB = 0
!        read (cline, *, iostat = input_status) &
!             cline_head, ibp_read, iunit1, iunit2,  &
!             imp1, imp2, imp1un, imp2un,             &
!             dist, factor, idummy, coef, ctype3, nHB
!        if(input_status > 0) then
!           error_message = 'read error =>' // cline
!           call util_error(ERROR%STOP_ALL, error_message)
!        end if
!
!        ibp = ibp + 1
!        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
!           imp1 = imp1 + ii
!           imp2 = imp2 + jj
!        else ! NATIVEINFO%ONE_BY_ONE
!           imp1 = imp1un + ii
!           imp2 = imp2un + jj
!        end if
!        irna_bp2unit(1, ibp) = imp2unit(imp1)
!        irna_bp2unit(2, ibp) = imp2unit(imp2)
!        irna_bp2mp(1, ibp) = imp1
!        irna_bp2mp(2, ibp) = imp2
!        rna_bp_nat(ibp) = dist
!        rna_bp_nat2(ibp) = dist**2
!        factor_rna_bp(ibp) = factor
!        irna_bp_dummy_mgo(ibp) = idummy
!        if (inmisc%flg_coef_from_ninfo) then
!           coef_rna_bp(ibp) = coef
!        endif
!
!        if (ctype3 /= '   ') then
!           if (nHB < 2) then
!              write(error_message,*) 'Error: in read_nativeinfo, invalid nHB',nHB
!              call util_error(ERROR%STOP_ALL, error_message)
!           else
!              nhb_bp(ibp) = nHB
!           endif
!        endif
!     end if
!     
!     ! ------------------------------------------------------------------
!     ! read basestack
!     if(cline(1:9) == 'basestack') then
!        ctype3 = '   '
!        read (cline, *, iostat = input_status) &
!             cline_head, ist_read, iunit1, iunit2,  &
!             imp1, imp2, imp1un, imp2un,             &
!             dist, factor, idummy, coef, ctype3
!        if(input_status > 0) then
!           error_message = 'read error =>' // cline
!           call util_error(ERROR%STOP_ALL, error_message)
!        end if
!
!        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
!           imp1 = imp1 + ii
!           imp2 = imp2 + jj
!        else ! NATIVEINFO%ONE_BY_ONE
!           imp1 = imp1un + ii
!           imp2 = imp2un + jj
!        end if
!
!        ist = ist + 1
!        irna_st2unit(1, ist) = imp2unit(imp1)
!        irna_st2unit(2, ist) = imp2unit(imp2)
!        irna_st2mp(1, ist) = imp1
!        irna_st2mp(2, ist) = imp2
!        rna_st_nat(ist) = dist
!        rna_st_nat2(ist) = dist**2
!        factor_rna_st(ist) = factor
!        irna_st_dummy_mgo(ist) = idummy
!        if (inmisc%flg_coef_from_ninfo) then
!           coef_rna_st(ist) = coef
!        endif
!     end if

     ! ------------------------------------------------------------------
     ! read basestack of DT model
     if(cline(1:7) == 'bs-dist') then
        ctype3 = '   '
        read (cline, *, iostat = input_status) &
             cline_head, ist_read, iunit1, iunit2,  &
             imp1, imp2, imp1un, imp2un,             &
             energy0, dist, coef, ctype3
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        if (ist_read < 1 .or. MXDTRNAST < ist_read) then
           error_message = 'ibs of bs-dist is out of range. You may need to increase MXDTRNAST in const_maxsize.F90'
           call util_error(ERROR%STOP_ALL, error_message)
        endif

        datacheck_dtrna_st_dihd(1,ist_read) = .False.
        datacheck_dtrna_st_dihd(2,ist_read) = .False.

        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + jj
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + jj
        end if

        if (imp1 > imp2) then
           imp_tmp = imp1
           imp1 = imp2
           imp2 = imp_tmp
        endif

        if (imp2 /= imp1+3) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        ibsdist = ibsdist + 1
        !!!! Caution: read this code carefully!
        !!!! The usage of ibsdist is unusual.
        !!!! Below, index is ist_read, but not ibsdist.

        idtrna_st2mp(1, ist_read) = imp1      ! B1
        idtrna_st2mp(2, ist_read) = imp2      ! B1
        idtrna_st2mp(3, ist_read) = imp1 - 2  ! P1
        idtrna_st2mp(4, ist_read) = imp1 - 1  ! S1
        idtrna_st2mp(5, ist_read) = imp1 + 1  ! P2
        idtrna_st2mp(6, ist_read) = imp1 + 2  ! S2
        idtrna_st2mp(7, ist_read) = imp1 + 4  ! P3

        dtrna_st_nat(1,ist_read) = dist
        if (inmisc%flg_coef_from_ninfo) then
           coef_dtrna_st(1, ist_read, :) = coef
           coef_dtrna_st(0, ist_read, :) = energy0
        endif

        ctype2(1:1) = ctype3(1:1)
        ctype2(2:2) = ctype3(3:3)
        idtrna_st2nn(ist_read) = ifunc_nn2id(ctype2)
     end if


     ! ------------------------------------------------------------------
     ! read hydrogen bond of DT model
     if(cline(1:7) == 'hb-dist') then
        if (inmisc%i_dtrna_model == 2015) then
           read (cline, *, iostat = input_status) &
                cline_head, ihb_read, iunit1, iunit2,  &
                imp1, imp2, imp1un, imp2un,             &
                energy0, dist, coef, ctype1, nHB
           select case (nHB)
           case (1)
              read (cline, *, iostat = input_status) &
                   cline_head, ihb_read, iunit1, iunit2,  &
                   imp1, imp2, imp1un, imp2un,             &
                   energy0, dist, coef, ctype1, nHB, a11, a12
           case (2)
              read (cline, *, iostat = input_status) &
                   cline_head, ihb_read, iunit1, iunit2,  &
                   imp1, imp2, imp1un, imp2un,             &
                   energy0, dist, coef, ctype1, nHB, a11, a12, a21, a22
           case (3)
              read (cline, *, iostat = input_status) &
                   cline_head, ihb_read, iunit1, iunit2,  &
                   imp1, imp2, imp1un, imp2un,             &
                   energy0, dist, coef, ctype1, nHB, a11, a12, a21, a22, a31, a32
           case default
              error_message = 'invalid nHB or atom names:' // cline
              call util_error(ERROR%STOP_ALL, error_message)
           endselect
        else
           read (cline, *, iostat = input_status) &
                cline_head, ihb_read, iunit1, iunit2,  &
                imp1, imp2, imp1un, imp2un,             &
                energy0, dist, coef
        endif
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        if (ihb_read < 1 .or. MXDTRNAHB < ihb_read) then
           error_message = 'ihb of hb-dist is out of range. You may need to increase MXDTRNAHB in const_maxsize.F90'
           call util_error(ERROR%STOP_ALL, error_message)
        endif

        datacheck_dtrna_hb_angl(1:2, ihb_read) = .False.
        datacheck_dtrna_hb_dihd(1:3, ihb_read) = .False.

        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + jj
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + jj
        end if

        ihbdist = ihbdist + 1
        if (imp1 < imp2) then
           idtrna_hb2mp(1, ihb_read) = imp1
           idtrna_hb2mp(2, ihb_read) = imp2
        else
           idtrna_hb2mp(1, ihb_read) = imp2
           idtrna_hb2mp(2, ihb_read) = imp1
        endif
        dtrna_hb_nat(1,ihb_read) = dist
        if (inmisc%flg_coef_from_ninfo) then
           coef_dtrna_hb(1, ihb_read) = coef
           coef_dtrna_hb(0, ihb_read) = energy0
        endif

        if (inmisc%i_dtrna_model == 2015) then
           if (imp2 < imp1) then
              atmp = a11
              a11 = a12
              a12 = atmp
              atmp = a21
              a21 = a22
              a22 = atmp
              atmp = a31
              a31 = a32
              a32 = atmp
           endif

           offset1 = imp2hbsite(1, imp1) - 1
           offset2 = imp2hbsite(1, imp2) - 1
           idtrna_hb2hbsite(1,1,ihb_read) = offset1 + site_local(cmp2seq(imp1), a11)
           idtrna_hb2hbsite(1,2,ihb_read) = offset2 + site_local(cmp2seq(imp2), a12)
           if (nHB > 1) then
              idtrna_hb2hbsite(2,1,ihb_read) = offset1 + site_local(cmp2seq(imp1), a21)
              idtrna_hb2hbsite(2,2,ihb_read) = offset2 + site_local(cmp2seq(imp2), a22)
           endif
           if (nHB > 2) then
              idtrna_hb2hbsite(3,1,ihb_read) = offset1 + site_local(cmp2seq(imp1), a31)
              idtrna_hb2hbsite(3,2,ihb_read) = offset2 + site_local(cmp2seq(imp2), a32)
           endif

           if (ctype1 == 'S') then
              flg_hb_tertiary(ihb_read) = .False.
           else if (ctype1 == 'T') then
              flg_hb_tertiary(ihb_read) = .True.
           else
              error_message = 'read error =>' // cline
              call util_error(ERROR%STOP_ALL, error_message)
           endif

        endif
     end if


     ! ------------------------------------------------------------------
     ! read basestack of DT model (2015)
     if(cline(1:8) == 'tbs-dist') then
        ctype3 = '   '
        read (cline, *, iostat = input_status) &
             cline_head, ist_read, iunit1, iunit2,  &
             imp1, imp2, imp1un, imp2un,             &
             energy0, dist, coef, ist_type1, ist_type2
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        if (ist_read < 1 .or. MXDTRNAHB < ist_read) then
           error_message = 'ist of st-dist is out of range. You may need to increase MXDTRNAST in const_maxsize.F90'
           call util_error(ERROR%STOP_ALL, error_message)
        endif

        datacheck_dtrna_tst_angl(1:2, ist_read) = .False.
        datacheck_dtrna_tst_dihd(1:3, ist_read) = .False.

        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + jj
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + jj
        end if

        if (imp1 > imp2) then
           imp_tmp = imp1
           imp1 = imp2
           imp2 = imp_tmp

           ist_type_tmp = ist_type1
           ist_type1 = ist_type2
           ist_type2 = ist_type_tmp
        endif

        itbsdist = itbsdist + 1
        !!!! Caution: read this code carefully!
        !!!! The usage of ibsdist is unusual.
        !!!! Below, index is ist_read, but not ibsdist.

        idtrna_tst2mp(1, ist_read) = imp1      ! B1
        idtrna_tst2mp(2, ist_read) = imp2      ! B2
        idtrna_tst2mp(3, ist_read) = imp1 - 1  ! S1
        idtrna_tst2mp(4, ist_read) = imp2 - 1  ! S2
        idtrna_tst2mp(5, ist_read) = imp1 + 1  ! P1(next)
        idtrna_tst2mp(6, ist_read) = imp2 + 1  ! P2(next)

        dtrna_tst_nat(1,ist_read) = dist
        if (inmisc%flg_coef_from_ninfo) then
           coef_dtrna_tst(1, ist_read) = coef
           coef_dtrna_tst(0, ist_read) = energy0
        endif

        if (ist_type1 < 0) then
           idtrna_tst2side(1,ist_read) = 1
        else
           idtrna_tst2side(1,ist_read) = 2
        endif
        if (ist_type2 < 0) then
           idtrna_tst2side(2,ist_read) = 1
        else
           idtrna_tst2side(2,ist_read) = 2
        endif

        if (abs(ist_type1) == 1) then
           flg_tst_exclusive(1,ist_read) = .False.
        else if (abs(ist_type1) == 2) then
           flg_tst_exclusive(1,ist_read) = .True.
        else
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        if (abs(ist_type2) == 1) then
           flg_tst_exclusive(2,ist_read) = .False.
        else if (abs(ist_type2) == 2) then
           flg_tst_exclusive(2,ist_read) = .True.
        else
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        endif

     end if

  end do

  ! ---------------------------------------------------------------------
  nbd = ibd 
  nfene = ifene
  nrouse = irouse
  nba = iba 
!  ndih = idih
  ncon = icon
  nLJ  = iLJ
  nwca  = iwca
  ncon_gauss  = icon_gauss
!  nrna_bp = ibp
!  nrna_st = ist
  ndtrna_st = ibsdist
  ndtrna_hb = ihbdist
  ndtrna_tst = itbsdist

  ! ---------------------------------------------------------------------
  ! check the input 
  do icon = 1, ncon
     if(icon2mp(1, icon) == 0 .or. icon2mp(2, icon) == 0) then
        error_message = 'Error: at contact in read_nativeinfo'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
  end do


  ! ------------------------------------------------------------------
  ! read again for bs-dihd, hb-angl, and hb-dihd
  ! ------------------------------------------------------------------
  rewind(lun)
  do
     read (lun, '(a256)', iostat = input_status) cline
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: cannot read intra coordinate in read_nativeinfo'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     ! ------------------------------------------------------------------
     ! read basestack dihedral of DT model
     if(cline(1:7) == 'bs-dihd') then
        ctype4 = '    '
        read (cline, *, iostat = input_status) &
             cline_head, ist_read, idih_read, iunit1, iunit2,  &
             imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
             dih, coef, ctype4
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        dih = dih * F_PI/180.0

        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + ii
           imp3 = imp3 + ii
           imp4 = imp4 + ii
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + ii
           imp3 = imp3un + ii
           imp4 = imp4un + ii
        end if

        if (ctype4 == 'PSPS') then
           if ((imp1 == idtrna_st2mp(3,ist_read) .and. &
                imp2 == idtrna_st2mp(4,ist_read) .and. &
                imp3 == idtrna_st2mp(5,ist_read) .and. &
                imp4 == idtrna_st2mp(6,ist_read)) .or. &
               (imp4 == idtrna_st2mp(3,ist_read) .and. &
                imp3 == idtrna_st2mp(4,ist_read) .and. &
                imp2 == idtrna_st2mp(5,ist_read) .and. &
                imp1 == idtrna_st2mp(6,ist_read)) ) then
              continue
           else
              error_message = 'read error =>' // cline
              call util_error(ERROR%STOP_ALL, error_message)
           endif
              
           dtrna_st_nat(2, ist_read) = dih
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_st(2, ist_read, :) = coef
           endif

           datacheck_dtrna_st_dihd(1,ist_read) = .True.

        else if (ctype4 == 'SPSP') then
           if ((imp1 == idtrna_st2mp(4,ist_read) .and. &
                imp2 == idtrna_st2mp(5,ist_read) .and. &
                imp3 == idtrna_st2mp(6,ist_read) .and. &
                imp4 == idtrna_st2mp(7,ist_read)) .or. &
               (imp4 == idtrna_st2mp(4,ist_read) .and. &
                imp3 == idtrna_st2mp(5,ist_read) .and. &
                imp2 == idtrna_st2mp(6,ist_read) .and. &
                imp1 == idtrna_st2mp(7,ist_read)) ) then
              continue
           else
              error_message = 'read error =>' // cline
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           dtrna_st_nat(3, ist_read) = dih
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_st(3, ist_read, :) = coef
           endif

           datacheck_dtrna_st_dihd(2,ist_read) = .True.

        else
           error_message = 'Error: unknown type of bs-dihd in read_nativeinfo'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
     end if

     ! ------------------------------------------------------------------
     ! read hydrogen-bond angle of DT model
     if(cline(1:7) == 'hb-angl') then
        read (cline, *, iostat = input_status) &
            cline_head, ihb_read, iba_read, iunit1, iunit2, &
            imp1, imp2, imp3, imp1un, imp2un, imp3un, &
            ba, coef
        if(input_status > 0) then
           write(*,*) input_status
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        ba = ba * F_PI/180.0

        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + ii
           imp3 = imp3 + ii
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + ii
           imp3 = imp3un + ii
        end if
        
        if (imp1 > imp3) then
           imp_tmp = imp1
           imp1 = imp3
           imp3 = imp_tmp
        endif

        if (imp2 == idtrna_hb2mp(1,ihb_read) .and. imp3 == idtrna_hb2mp(2,ihb_read)) then
           if (imp1 == idtrna_hb2mp(3,ihb_read)) then
              continue
           else if (idtrna_hb2mp(3,ihb_read) == 0) then
              idtrna_hb2mp(3,ihb_read) = imp1
           else
              error_message = 'read error =>' // cline
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           dtrna_hb_nat(2,ihb_read) = ba
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_hb(2,ihb_read) = coef
           endif
           
           datacheck_dtrna_hb_angl(1, ihb_read) = .True.

        else if (imp1 == idtrna_hb2mp(1,ihb_read) .and. imp2 == idtrna_hb2mp(2,ihb_read)) then
           if (imp3 == idtrna_hb2mp(4,ihb_read)) then
              continue
           else if (idtrna_hb2mp(4,ihb_read) == 0) then
              idtrna_hb2mp(4,ihb_read) = imp3
           else
              error_message = 'read error =>' // cline
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           dtrna_hb_nat(3,ihb_read) = ba
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_hb(3,ihb_read) = coef
           endif

           datacheck_dtrna_hb_angl(2, ihb_read) = .True.

        else
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        endif
     end if

     ! ------------------------------------------------------------------
     ! read hydrogen-bond dihedral of DT model
     if(cline(1:7) == 'hb-dihd') then
        read (cline, *, iostat = input_status) &
             cline_head, ihb_read, idih_read, iunit1, iunit2,  &
             imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
             dih, coef
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        dih = dih * F_PI/180.0

        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + ii
           imp3 = imp3 + ii
           imp4 = imp4 + ii
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + ii
           imp3 = imp3un + ii
           imp4 = imp4un + ii
        end if

        if (imp1 > imp4) then
           imp_tmp = imp1
           imp1 = imp4
           imp4 = imp_tmp
           imp_tmp = imp2
           imp2 = imp3
           imp3 = imp_tmp
        endif

        if (imp2 == idtrna_hb2mp(1,ihb_read) .and. imp3 == idtrna_hb2mp(2,ihb_read)) then
           if (imp1 == idtrna_hb2mp(3, ihb_read)) then
              continue
           else if (idtrna_hb2mp(3, ihb_read) == 0) then
              idtrna_hb2mp(3,ihb_read) = imp1
           else
              error_message = 'read error =>' // cline
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           if (imp4 == idtrna_hb2mp(4, ihb_read)) then
              continue
           else if (idtrna_hb2mp(4, ihb_read) == 0) then
              idtrna_hb2mp(4,ihb_read) = imp4
           else
              error_message = 'read error =>' // cline
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           dtrna_hb_nat(4, ihb_read) = dih
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_hb(4, ihb_read) = coef
           endif

           datacheck_dtrna_hb_dihd(1, ihb_read) = .True.

        else if (imp3 == idtrna_hb2mp(1,ihb_read) .and. imp4 == idtrna_hb2mp(2,ihb_read)) then
           if (imp2 == idtrna_hb2mp(3, ihb_read)) then
              continue
           else if (idtrna_hb2mp(3, ihb_read) == 0) then
              idtrna_hb2mp(3,ihb_read) = imp2
           else
              error_message = 'read error =>' // cline
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           idtrna_hb2mp(5,ihb_read) = imp1
           dtrna_hb_nat(5, ihb_read) = dih
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_hb(5, ihb_read) = coef
           endif

           datacheck_dtrna_hb_dihd(2, ihb_read) = .True.

        else if (imp1 == idtrna_hb2mp(1,ihb_read) .and. imp2 == idtrna_hb2mp(2,ihb_read)) then
           if (imp3 == idtrna_hb2mp(4,ihb_read)) then
              continue
           else if (idtrna_hb2mp(4,ihb_read) == 0) then
              idtrna_hb2mp(4,ihb_read) = imp3
           else
              error_message = 'read error =>' // cline
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           idtrna_hb2mp(6,ihb_read) = imp4
           dtrna_hb_nat(6, ihb_read) = dih
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_hb(6, ihb_read) = coef
           endif

           datacheck_dtrna_hb_dihd(3, ihb_read) = .True.

        else
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        endif
     endif

     ! ------------------------------------------------------------------
     ! read nonlocal stacking angle of DT model (2015)
     if(cline(1:8) == 'tbs-angl') then
        read (cline, *, iostat = input_status) &
            cline_head, ist_read, iba_read, iunit1, iunit2, &
            imp1, imp2, imp3, imp1un, imp2un, imp3un, &
            ba, coef
        if(input_status > 0) then
           write(*,*) input_status
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        ba = ba * F_PI/180.0

        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + ii
           imp3 = imp3 + ii
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + ii
           imp3 = imp3un + ii
        end if
        
        if (imp1 > imp3) then
           imp_tmp = imp1
           imp1 = imp3
           imp3 = imp_tmp
        endif

        if (imp1 == idtrna_tst2mp(3,ist_read) .and.& 
            imp2 == idtrna_tst2mp(1,ist_read) .and.&
            imp3 == idtrna_tst2mp(2,ist_read) ) then
           dtrna_tst_nat(2,ist_read) = ba
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_tst(2,ist_read) = coef
           endif

           datacheck_dtrna_tst_angl(1, ist_read) = .True.

        else if (imp1 == idtrna_tst2mp(1,ist_read) .and.&
                 imp2 == idtrna_tst2mp(2,ist_read) .and.&
                 imp3 == idtrna_tst2mp(4,ist_read) ) then
           dtrna_tst_nat(3,ist_read) = ba
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_tst(3,ist_read) = coef
           endif

           datacheck_dtrna_tst_angl(2, ist_read) = .True.

        else
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        endif
     end if

     ! ------------------------------------------------------------------
     ! read nonlocal stacking dihedral of DT model (2015)
     if(cline(1:8) == 'tbs-dihd') then
        read (cline, *, iostat = input_status) &
             cline_head, ist_read, idih_read, iunit1, iunit2,  &
             imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
             dih, coef
        if(input_status > 0) then
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        dih = dih * F_PI/180.0

        if(i_ninfo_type == NATIVEINFO%ALL_IN_ONE) then
           imp1 = imp1 + ii
           imp2 = imp2 + ii
           imp3 = imp3 + ii
           imp4 = imp4 + ii
        else ! NATIVEINFO%ONE_BY_ONE
           imp1 = imp1un + ii
           imp2 = imp2un + ii
           imp3 = imp3un + ii
           imp4 = imp4un + ii
        end if

        if (imp1 > imp4) then
           imp_tmp = imp1
           imp1 = imp4
           imp4 = imp_tmp
           imp_tmp = imp2
           imp2 = imp3
           imp3 = imp_tmp
        endif

        if (imp1 == idtrna_tst2mp(3,ist_read) .and.&
            imp2 == idtrna_tst2mp(1,ist_read) .and.&
            imp3 == idtrna_tst2mp(2,ist_read) .and.&
            imp4 == idtrna_tst2mp(4,ist_read) ) then
           dtrna_tst_nat(4, ist_read) = dih
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_tst(4, ist_read) = coef
           endif

           datacheck_dtrna_tst_dihd(1, ist_read) = .True.

        elseif (imp1 == idtrna_tst2mp(5,ist_read) .and.&
                imp2 == idtrna_tst2mp(3,ist_read) .and.&
                imp3 == idtrna_tst2mp(1,ist_read) .and.&
                imp4 == idtrna_tst2mp(2,ist_read) ) then
           dtrna_tst_nat(5,ist_read) = dih
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_tst(5, ist_read) = coef
           endif

           datacheck_dtrna_tst_dihd(2, ist_read) = .True.

        elseif (imp1 == idtrna_tst2mp(1,ist_read) .and.&
                imp2 == idtrna_tst2mp(2,ist_read) .and.&
                imp3 == idtrna_tst2mp(4,ist_read) .and.&
                imp4 == idtrna_tst2mp(6,ist_read) ) then
           dtrna_tst_nat(6,ist_read) = dih
           if (inmisc%flg_coef_from_ninfo) then
              coef_dtrna_tst(6, ist_read) = coef
           endif

           datacheck_dtrna_tst_dihd(3, ist_read) = .True.

        else
           error_message = 'read error =>' // cline
           call util_error(ERROR%STOP_ALL, error_message)
        endif
     endif
  enddo

  !!! Check
  do ihb_read = 1, ndtrna_hb
     if (datacheck_dtrna_hb_angl(1,ihb_read) .and. &
         datacheck_dtrna_hb_angl(2,ihb_read) ) then
        continue
     else
        write(error_message,*) 'Error in ninfo: following ihb does not have angle: ', ihb_read
        call util_error(ERROR%STOP_ALL, error_message)
     endif
     if (datacheck_dtrna_hb_dihd(1,ihb_read) .and. &
         datacheck_dtrna_hb_dihd(2,ihb_read) .and. &
         datacheck_dtrna_hb_dihd(3,ihb_read)  ) then
        continue
     else
        write(error_message,*) 'Error in ninfo: following ihb does not have dihedral: ', ihb_read
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  enddo
  do ist_read = 1, ndtrna_st
     if (datacheck_dtrna_st_dihd(1,ist_read) .and. &
         datacheck_dtrna_st_dihd(2,ist_read) ) then
        continue
     else
        write(error_message,*) 'Error in ninfo: following ist does not have dihedral: ', ist_read
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  enddo
  do ist_read = 1, ndtrna_tst
     if (datacheck_dtrna_tst_angl(1,ist_read) .and. &
         datacheck_dtrna_tst_angl(2,ist_read) ) then
        continue
     else
        write(error_message,*) 'Error in ninfo: following itst does not have angle: ', ist_read
        call util_error(ERROR%STOP_ALL, error_message)
     endif
     if (datacheck_dtrna_tst_dihd(1,ist_read) .and. &
         datacheck_dtrna_tst_dihd(2,ist_read) .and. &
         datacheck_dtrna_tst_dihd(3,ist_read) ) then
        continue
     else
        write(error_message,*) 'Error in ninfo: following itst does not have dihedral: ', ist_read
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  enddo
   

  ! ---------------------------------------------------------------------


!############################################################################
contains
   integer function str2bondtype(c2)
      use const_index
      implicit none
      character(2), intent(in) :: c2
      if (c2 == 'pp') then
         str2bondtype = BDTYPE%PRO
      else if (c2 == 'PS') then
         str2bondtype = BDTYPE%RNA_PS
      else if (c2 == 'SP') then
         str2bondtype = BDTYPE%RNA_SP
      else if (c2 == 'SA' .OR. c2 == 'SG' .OR. c2 == 'SR') then
         str2bondtype = BDTYPE%RNA_SR
      else if (c2 == 'SU' .OR. c2 == 'SC' .OR. c2 == 'SY') then
         str2bondtype = BDTYPE%RNA_SY
      else
         str2bondtype = BDTYPE%VOID
         error_message = 'Error: in read_nativeinfo, unknown bondtype'//c2
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   endfunction str2bondtype

   integer function str2angletype(c3)
      use const_index
      implicit none
      character(3), intent(in) :: c3
      if (c3 == 'ppp') then
         str2angletype = BATYPE%PRO
      else if (c3 == 'PSP') then
         str2angletype = BATYPE%RNA_PSP
      else if (c3 == 'SPS') then
         str2angletype = BATYPE%RNA_SPS
      else if (c3 == 'ASP' .OR. c3 == 'GSP' .OR. c3 == 'RSP') then
         str2angletype = BATYPE%RNA_RSP
      else if (c3 == 'USP' .OR. c3 == 'CSP' .OR. c3 == 'YSP') then
         str2angletype = BATYPE%RNA_YSP
      else if (c3 == 'PSA' .OR. c3 == 'PSG' .OR. c3 == 'PSR') then
         str2angletype = BATYPE%RNA_PSR
      else if (c3 == 'PSU' .OR. c3 == 'PSC' .OR. c3 == 'PSY') then
         str2angletype = BATYPE%RNA_PSY
      else
         str2angletype = BATYPE%VOID
         error_message = 'Error: in read_nativeinfo, unknown angletype'//c3
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   endfunction str2angletype

   integer function str2dihtype(c4)
      use const_index
      implicit none
      character(4), intent(in) :: c4
      if (c4 == 'pppp') then
         str2dihtype = DIHTYPE%PRO
      else if (c4 == 'PSPS') then
         str2dihtype = DIHTYPE%RNA_PSPS
      else if (c4 == 'SPSP') then
         str2dihtype = DIHTYPE%RNA_SPSP
      else if (c4 == 'ASPS' .OR. c4 == 'GSPS' .OR. c4 == 'RSPS') then
         str2dihtype = DIHTYPE%RNA_RSPS
      else if (c4 == 'USPS' .OR. c4 == 'CSPS' .OR. c4 == 'YSPS') then
         str2dihtype = DIHTYPE%RNA_YSPS
      else if (c4 == 'SPSA' .OR. c4 == 'SPSG' .OR. c4 == 'SPSR') then
         str2dihtype = DIHTYPE%RNA_SPSR
      else if (c4 == 'SPSU' .OR. c4 == 'SPSC' .OR. c4 == 'SPSY') then
         str2dihtype = DIHTYPE%RNA_SPSY
      else if (c4 == 'ASSA' .OR. c4 == 'ASSG' .OR. c4 == 'ASSR' .OR. &
               c4 == 'GSSA' .OR. c4 == 'GSSG' .OR. c4 == 'GSSR' .OR. &
               c4 == 'RSSA' .OR. c4 == 'RSSG' .OR. c4 == 'RSSR' .OR. &
               c4 == 'ASSU' .OR. c4 == 'ASSC' .OR. c4 == 'ASSY' .OR. &
               c4 == 'GSSU' .OR. c4 == 'GSSC' .OR. c4 == 'GSSY' .OR. &
               c4 == 'RSSU' .OR. c4 == 'RSSC' .OR. c4 == 'RSSY' .OR. &
               c4 == 'USSA' .OR. c4 == 'USSG' .OR. c4 == 'USSR' .OR. &
               c4 == 'CSSA' .OR. c4 == 'CSSG' .OR. c4 == 'CSSR' .OR. &
               c4 == 'YSSA' .OR. c4 == 'YSSG' .OR. c4 == 'YSSR' .OR. &
               c4 == 'USSU' .OR. c4 == 'USSC' .OR. c4 == 'USSY' .OR. &
               c4 == 'CSSU' .OR. c4 == 'CSSC' .OR. c4 == 'CSSY' .OR. &
               c4 == 'YSSU' .OR. c4 == 'YSSC' .OR. c4 == 'YSSY' ) then
         ! These are "stack dihedral" supported in previous version.
         str2dihtype = DIHTYPE%VOID
      else
         str2dihtype = DIHTYPE%VOID
         error_message = 'Error: in read_nativeinfo, unknown dihtype'//c4
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   endfunction str2dihtype

   integer function str2gotype(c3)
      use const_index
      implicit none
      character(3), intent(in) :: c3

      if (c3 == 'p-p') then
         str2gotype = CONTYPE%PRO_PRO
      else if (c3 == 'p-P') then
         str2gotype = CONTYPE%PRO_RP
      else if (c3 == 'p-S') then
         str2gotype = CONTYPE%PRO_RS
      else if (c3 == 'p-A' .OR. c3 == 'p-G' .OR. c3 == 'p-R' .OR. &
               c3 == 'p-U' .OR. c3 == 'p-C' .OR. c3 == 'p-Y' ) then
         str2gotype = CONTYPE%PRO_RB
      else if (c3 == 'P-P') then
         str2gotype = CONTYPE%RP_RP
      else if (c3 == 'P-S') then
         str2gotype = CONTYPE%RP_RS
      else if (c3 == 'P-A' .OR. c3 == 'P-G' .OR. c3 == 'P-R' .OR. &
               c3 == 'P-U' .OR. c3 == 'P-C' .OR. c3 == 'P-Y' ) then
         str2gotype = CONTYPE%RP_RB
      else if (c3 == 'S-S') then
         str2gotype = CONTYPE%RS_RS
      else if (c3 == 'S-A' .OR. c3 == 'S-G' .OR. c3 == 'S-R' .OR. &
               c3 == 'S-U' .OR. c3 == 'S-C' .OR. c3 == 'S-Y' ) then
         str2gotype = CONTYPE%RS_RB
      else if (c3 == 'A-A' .OR. c3 == 'A-G' .OR. c3 == 'A-R' .OR. &
               c3 == 'A-U' .OR. c3 == 'A-C' .OR. c3 == 'A-Y' .OR. &
               c3 == 'U-G' .OR. c3 == 'U-R' .OR. &
               c3 == 'U-U' .OR. c3 == 'U-C' .OR. c3 == 'U-Y' .OR. &
               c3 == 'G-G' .OR. c3 == 'G-R' .OR. &
               c3 == 'G-C' .OR. c3 == 'G-Y' .OR. &
               c3 == 'C-R' .OR. c3 == 'C-C' .OR. c3 == 'C-Y' .OR. &
               c3 == 'R-R' .OR. c3 == 'R-Y' .OR. c3 == 'Y-Y' ) then
         str2gotype = CONTYPE%RB_RB
      else
         error_message = 'Error: in read_nativeinfo, unknown contact type'//c3
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   endfunction str2gotype

   integer function site_local(res, atom)
      implicit none
      character(3), intent(in) :: res
      character(4), intent(in) :: atom
      
      site_local = 0
      if (atom(1:3) == "O2'") then
         site_local = 1
      else if (atom(1:3) == "O4'") then
         site_local = 2
      else if (atom(1:3) == "O3'") then
         site_local = 3
      else if (atom(1:3) == "O5'") then
         site_local = 4
      else if (atom(1:3) == "OP1") then
         site_local = 1
      else if (atom(1:3) == "OP2") then
         site_local = 2
      else if (atom(1:2) == "N1") then
         site_local = 1
      else if (atom(1:2) == "O2") then
         site_local = 1
      else if (atom(1:2) == "N2") then
         site_local = 2
      else if (atom(1:2) == "N3") then
         if (res(1:2) == "RG") then
            site_local = 3
         else if (res(1:2) == "RA" .or. res(1:2) == "RC" .or. res(1:2) == "RU") then
            site_local = 2
         endif
      else if (atom(1:2) == "O4") then
         site_local = 3
      else if (atom(1:2) == "N4") then
         site_local = 3
      else if (atom(1:2) == "N6") then
         site_local = 3
      else if (atom(1:2) == "O6") then
         site_local = 4
      else if (atom(1:2) == "N7") then
         if (res(1:2) == "RA") then
            site_local = 4
         else if (res(1:2) == "RG") then
            site_local = 5
         endif
      endif

      if (site_local == 0) then
         error_message = 'Error: in hbsite_local (read_nativeinfo): '//res//','//atom
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   endfunction site_local 

end subroutine read_nativeinfo
