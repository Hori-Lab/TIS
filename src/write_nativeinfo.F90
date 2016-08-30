! write_nativeinfo
!> @brief Write nativeinfo to *.ninfo file

subroutine write_nativeinfo(lunout)
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inpro
  use var_struct, only : nunit_all, imp2unit, lunit2mp, &
                         nbd, ibd2mp, bd_nat, cmp2atom, &
                         factor_bd, coef_bd, correct_bd_mgo,  &
                         nfene, ifene2mp, fene_nat, coef_fene, dist2_fene,  &
                         nba, iba2mp, ba_nat, factor_ba, coef_ba, correct_ba_mgo,&
!                         ndih, idih2mp, dih_nat,    &
!                         factor_dih, coef_dih, correct_dih_mgo, &
                         ncon, icon2mp, factor_go,       &
                         coef_go, icon_dummy_mgo, go_nat, ncon_unit, &
                         nLJ, iLJ2mp, coef_LJ, LJ_nat,  &
                         iclass_unit, ibd2type, iba2type, idih2type, icon2type, &
!                         nrna_bp,nrna_bp_unit, nrna_st, rna_bp_nat, &
!                         irna_bp2mp, coef_rna_bp, factor_rna_bp, nhb_bp, &
!                         irna_st2mp, coef_rna_st, factor_rna_st, rna_st_nat, &
!                         coef_aicg13_gauss, wid_aicg13_gauss, aicg13_nat, factor_aicg13, & ! AICG2
!                         coef_aicg14_gauss, wid_aicg14_gauss, aicg14_nat, factor_aicg14, & ! AICG2
!                         coef_dih_gauss, wid_dih_gauss, & ! AICG2
                         ndtrna_st, idtrna_st2mp, dtrna_st_nat, coef_dtrna_st, &
                         ndtrna_hb, idtrna_hb2mp, dtrna_hb_nat, coef_dtrna_hb
  use mpiconst

  implicit none

  integer, intent(in) :: lunout

  integer :: iunit, junit
  integer :: imp1, imp2, imp3, imp4, iunit1, iunit2
  integer :: imp1un, imp2un, imp3un, imp4un
  integer :: ibd, iba, idih, icon!, ibp
  real(PREC) :: dfcontact
  character(CARRAY_MSG_ERROR) :: error_message
  integer, parameter :: IREP = 1

  ! ------------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  ! -------------------------------------------------------------------
  ! write the bond length
  if (nbd > 0) then
     write (lunout, '(a)') '<<<< native bond length '
     write (lunout, '(a)') '** coef_bd(kcal/mol) = factor_bd * correct_bd_mgo * cbd * energy_unit_protein'
     write (lunout, '(a)', ADVANCE='NO') '**      ibd iunit1-iunit2   imp1 - imp2 imp1un-imp2un'
     write (lunout, '(a)')               '      bd_nat    factor_bd  correct_mgo      coef_bd'
       
     do ibd = 1, nbd
        imp1 = ibd2mp(1, ibd)
        imp2 = ibd2mp(2, ibd)
        iunit1 = imp2unit(imp1)
        iunit2 = iunit1
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit1) + 1
        write (lunout, "(a4, 7(1xi6), 4(1xf12.4))", ADVANCE='NO') &
             'bond', ibd, iunit1, iunit2, imp1, imp2, imp1un, imp2un, &
             bd_nat(ibd), factor_bd(ibd), correct_bd_mgo(ibd), coef_bd(1, ibd)
        if (iclass_unit(iunit1) == CLASS%PRO) then
           write(lunout, '(a3)') ' pp'
        else if (iclass_unit(iunit1) == CLASS%RNA) then
           write(lunout, '(a3)') bondtype2str()
        else
           write(lunout, '(a)') ''
        endif
     end do
        
     write (lunout, '(a4)') '>>>>'
     write (lunout, '(a)') ''
  endif  ! nbd > 0

  
  ! -------------------------------------------------------------------
  ! FENE
  if (nfene > 0) then
     write (lunout, '(a)') '<<<< FENE'
     write (lunout, '(a)', ADVANCE='NO') '**      ibd iunit1-iunit2   imp1 - imp2 imp1un-imp2un'
     write (lunout, '(a)')               '        dist        dist2         coef'
       
     do ibd = 1, nfene
        imp1 = ifene2mp(1, ibd)
        imp2 = ifene2mp(2, ibd)
        iunit1 = imp2unit(imp1)
        iunit2 = iunit1
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit1) + 1
        write (lunout, "(a4, 7(1xi6), 3(1xf12.4))") &
             'FENE', ibd, iunit1, iunit2, imp1, imp2, imp1un, imp2un, &
             fene_nat(ibd), dist2_fene(ibd), coef_fene(ibd)
     end do
        
     write (lunout, '(a4)') '>>>>'
     write (lunout, '(a)') ''
  endif  ! nfene > 0

  
  ! -------------------------------------------------------------------
  ! write the bond angle
  if (nba > 0) then
     write (lunout, '(a)') '<<<< native bond angles '
     write (lunout, '(a)') '** coef_ba(kcal/mol) = factor_ba * correct_ba_mgo * cba * energy_unit_protein'
     write (lunout, '(a)', ADVANCE='NO') '**      iba iunit1-iunit2   imp1 - imp2 - imp3'
     write (lunout, '(a)')               ' imp1un-imp2un-imp3un      ba_nat    factor_ba  correct_mgo      coef_ba'
     do iba = 1, nba
        imp1 = iba2mp(1, iba)
        imp2 = iba2mp(2, iba)
        imp3 = iba2mp(3, iba)
        iunit1 = imp2unit(imp1)
        iunit2 = iunit1
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit1) + 1
        imp3un = imp3 - lunit2mp(1, iunit1) + 1
        write (lunout, "(a4, 9(1xi6), 4(1xf12.4))", ADVANCE='NO') &
             'angl', iba, iunit1, iunit2, &
             imp1, imp2, imp3, imp1un, imp2un, imp3un, &
             ba_nat(iba) * 180.0e0_PREC / F_PI, &
             factor_ba(iba), correct_ba_mgo(iba), coef_ba(1, iba)
        if (iclass_unit(iunit1) == CLASS%PRO) then
           write(lunout, '(a4)') ' ppp'
        else if (iclass_unit(iunit1) == CLASS%RNA) then
           write(lunout, '(a4)') angletype2str()
        else 
           write(lunout, '(a)') ''
        endif
     end do
     write (lunout, '(a4)') '>>>>'
     write (lunout, '(a)') ''
  endif  ! nba > 0

!  ! -------------------------------------------------------------------
!  if (inmisc%force_flag_local(LINTERACT%L_AICG2) .or. &
!      inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
!  ! write the aicg13
!  if (nba > 0) then
!     write (lunout, '(a)') '<<<< 1-3 contacts with L_AICG2 or L_AICG2_PLUS'
!     write (lunout, '(a)') '** coef_aicg13_gauss(kcal/mol) = factor_aicg13 * correct_ba_mgo * coef_aicg13_gauss * energy_unit_protein'
!     write (lunout, '(a)', ADVANCE='NO') '**      iba iunit1-iunit2   imp1 - imp2 - imp3'
!     write (lunout, '(a)')               ' imp1un-imp2un-imp3un  aicg13_nat  factor_aicg13  correct_mgo  coef_aicg13_gauss wid_aicg13_gauss'
!     do iba = 1, nba
!        imp1 = iba2mp(1, iba)
!        imp2 = iba2mp(2, iba)
!        imp3 = iba2mp(3, iba)
!        iunit1 = imp2unit(imp1)
!        iunit2 = iunit1
!        imp1un = imp1 - lunit2mp(1, iunit1) + 1
!        imp2un = imp2 - lunit2mp(1, iunit1) + 1
!        imp3un = imp3 - lunit2mp(1, iunit1) + 1
!        write (lunout, "(a6, 9(1xi6), 5(1xf12.4))", ADVANCE='NO') &
!             'aicg13', iba, iunit1, iunit2, &
!             imp1, imp2, imp3, imp1un, imp2un, imp3un, &
!             aicg13_nat(iba), &
!             factor_aicg13(iba), correct_ba_mgo(iba), coef_aicg13_gauss(iba), wid_aicg13_gauss(iba)
!        if (iclass_unit(iunit1) == CLASS%PRO) then
!           write(lunout, '(a4)') ' ppp'
!        else if (iclass_unit(iunit1) == CLASS%RNA) then
!           write(lunout, '(a4)') angletype2str()
!        else
!           write(lunout, '(a)') ''
!        endif
!     end do
!     write (lunout, '(a4)') '>>>>'
!     write (lunout, '(a)') ''
!  endif  ! nba > 0
!  endif
  
!  ! -------------------------------------------------------------------
!  ! write the dihedral angle
!  if (ndih > 0) then
!     write (lunout, '(a)') '<<<< native dihedral angles '
!     write (lunout, '(a)') '** coef_dih1(kcal/mol) = factor_dih * correct_dih_mgo * cdih_1 * energy_unit_protein'
!     write (lunout, '(a)') '** coef_dih3(kcal/mol) = factor_dih * correct_dih_mgo * cdih_3 * energy_unit_protein'
!     write (lunout, '(a)', ADVANCE='NO') '**     idih iunit1-iunit2   imp1 - imp2 - imp3 - imp4'
!     write (lunout, '(a)') ' imp1un-imp2un-imp3un-imp4un      dih_nat   factor_dih  correct_mgo   coef_dih_1   coef_dih_3'
!   
!     do idih = 1, ndih   
!        imp1 = idih2mp(1, idih)
!        imp2 = idih2mp(2, idih)
!        imp3 = idih2mp(3, idih)
!        imp4 = idih2mp(4, idih)
!        iunit1 = imp2unit(imp1)
!        iunit2 = iunit1
!        imp1un = imp1 - lunit2mp(1, iunit1) + 1
!        imp2un = imp2 - lunit2mp(1, iunit1) + 1
!        imp3un = imp3 - lunit2mp(1, iunit1) + 1
!        imp4un = imp4 - lunit2mp(1, iunit1) + 1
!        write (lunout, "(a4, 11(1xi6), 5(1xf12.4))", ADVANCE='NO') &
!             'dihd', idih, iunit1, iunit2, imp1, imp2, imp3, imp4, &
!             imp1un, imp2un, imp3un, imp4un, &
!             dih_nat(idih) * 180.0e0_PREC / F_PI, &
!             factor_dih(idih), correct_dih_mgo(idih), &
!             coef_dih(1, idih), coef_dih(2, idih)
!        if (iclass_unit(iunit1) == CLASS%PRO) then
!           write(lunout, '(a5)') ' pppp'
!        else if (iclass_unit(iunit1) == CLASS%RNA) then
!           write(lunout, '(a5)') dihtype2str()
!        else
!           write(lunout, '(a)') ''
!        endif
!     end do
!     write (lunout, '(a4)') '>>>>'
!     write (lunout, '(a)') ''
!  endif  ! ndih > 0
  
!  ! -------------------------------------------------------------------
!  if (inmisc%force_flag_local(LINTERACT%L_AICG2)) then
!  ! write the aicg14
!  if (ndih > 0) then
!     write (lunout, '(a)') '<<<< 1-4 contacts with L_AICG2 '
!     write (lunout, '(a)') '** coef_aicg14_gauss(kcal/mol) = factor_aicg14 * correct_dih_mgo * coef_aicg14(kcal/mol) * energy_unit_protein'
!     write (lunout, '(a)', ADVANCE='NO') '**     idih iunit1-iunit2   imp1 - imp2 - imp3 - imp4'
!     write (lunout, '(a)') ' imp1un-imp2un-imp3un-imp4un   aicg14_nat factor_aicg14  correct_mgo  coef_aicg14_gauss  wid_aicg14_gauss'
!
!     do idih = 1, ndih
!        imp1 = idih2mp(1, idih)
!        imp2 = idih2mp(2, idih)
!        imp3 = idih2mp(3, idih)
!        imp4 = idih2mp(4, idih)
!        iunit1 = imp2unit(imp1)
!        iunit2 = iunit1
!        imp1un = imp1 - lunit2mp(1, iunit1) + 1
!        imp2un = imp2 - lunit2mp(1, iunit1) + 1
!        imp3un = imp3 - lunit2mp(1, iunit1) + 1
!        imp4un = imp4 - lunit2mp(1, iunit1) + 1
!        write (lunout, "(a6, 11(1xi6), 5(1xf12.4))", ADVANCE='NO') &
!             'aicg14', idih, iunit1, iunit2, imp1, imp2, imp3, imp4, &
!             imp1un, imp2un, imp3un, imp4un, &
!             aicg14_nat(idih), &
!             factor_aicg14(idih), correct_dih_mgo(idih), &
!             coef_aicg14_gauss(idih), wid_aicg14_gauss(idih)
!        if (iclass_unit(iunit1) == CLASS%PRO) then
!           write(lunout, '(a5)') ' pppp'
!        else if (iclass_unit(iunit1) == CLASS%RNA) then
!           write(lunout, '(a5)') dihtype2str()
!        else
!           write(lunout, '(a)') ''
!        endif
!     end do
!     write (lunout, '(a4)') '>>>>'
!     write (lunout, '(a)') ''
!  endif  ! ndih > 0
!  endif
!
!  ! -------------------------------------------------------------------
!  if (inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then
!  ! write the aicg14
!  if (ndih > 0) then
!     write (lunout, '(a)') '<<<< <<<< 1-4 contacts with L_AICG2_PLUS'
!     write (lunout, '(a)') '** coef_dih_gauss(kcal/mol) = factor_aicg14 * correct_dih_mgo * coef_aicg14(kcal/mol) * energy_unit_protein'
!     write (lunout, '(a)', ADVANCE='NO') '**     idih iunit1-iunit2   imp1 - imp2 - imp3 - imp4'
!     write (lunout, '(a)') ' imp1un-imp2un-imp3un-imp4un   dih_nat factor_aicg14  correct_mgo  coef_dih_gauss  wid_dih_gauss'
!
!     do idih = 1, ndih
!        imp1 = idih2mp(1, idih)
!        imp2 = idih2mp(2, idih)
!        imp3 = idih2mp(3, idih)
!        imp4 = idih2mp(4, idih)
!        iunit1 = imp2unit(imp1)
!        iunit2 = iunit1
!        imp1un = imp1 - lunit2mp(1, iunit1) + 1
!        imp2un = imp2 - lunit2mp(1, iunit1) + 1
!        imp3un = imp3 - lunit2mp(1, iunit1) + 1
!        imp4un = imp4 - lunit2mp(1, iunit1) + 1
!        write (lunout, "(a7, 11(1xi6), 5(1xf12.4))", ADVANCE='NO') &
!             'aicgdih', idih, iunit1, iunit2, imp1, imp2, imp3, imp4, &
!             imp1un, imp2un, imp3un, imp4un, &
!             dih_nat(idih) * 180.0e0_PREC / F_PI, &
!             factor_aicg14(idih), correct_dih_mgo(idih), &
!             coef_dih_gauss(idih), wid_dih_gauss(idih)
!        if (iclass_unit(iunit1) == CLASS%PRO) then
!           write(lunout, '(a5)') ' pppp'
!        else if (iclass_unit(iunit1) == CLASS%RNA) then
!           write(lunout, '(a5)') dihtype2str()
!        else
!           write(lunout, '(a)') ''
!        endif
!     end do
!     write (lunout, '(a4)') '>>>>'
!     write (lunout, '(a)') ''
!  endif  ! ndih > 0
!  endif

  ! ------------------------------------------------------------------
  ! write the go interaction
  if (ncon > 0) then
     dfcontact = inpro%dfcontact
!     if(inenm%i_enm == 1) dfcontact = inenm%dfcontact_enm

     write (lunout, '(a)') '<<<< native contact '
     write (lunout, '(a, i6)') '** total_contact = ', ncon
     write (lunout, '(a, f10.2, a)') '** definition_of_contact = ', dfcontact, ' A'
     write (lunout, '(a)') '** coef_go(kcal/mol) = factor_go * icon_dummy_mgo * cgo1210 * energy_unit_protein'
     write (lunout, '(a)') ''
     
     do iunit = 1, nunit_all
        do junit = iunit, nunit_all
           write (lunout, '(a, i6, a, i6)') '** contact between unit ', iunit, ' and ', junit
           !write (lunout, '(a, i6)') '** total_contact_unit = ', ncon_unit(iunit, junit) - nrna_bp_unit(iunit,junit)
           write (lunout, '(a, i6)') '** total_contact_unit = ', ncon_unit(iunit, junit) 
           write (lunout, '(a)', ADVANCE='NO') '**        icon iunit1-iunit2   imp1 - imp2 imp1un-imp2un'
           write (lunout, '(a)') '      go_nat   factor_go  dummy     coef_go'
           do icon = 1, ncon
              imp1 = icon2mp(1, icon)
              imp2 = icon2mp(2, icon)
              iunit1 = imp2unit(imp1)
              iunit2 = imp2unit(imp2)
              imp1un = imp1 - lunit2mp(1, iunit1) + 1
              imp2un = imp2 - lunit2mp(1, iunit2) + 1
                 
              if(iunit == iunit1 .and. junit == iunit2) then
                 if (iclass_unit(iunit) == CLASS%RNA .OR. iclass_unit(junit) == CLASS%RNA) then
                    if (icon2type(icon) == CONTYPE%RNA_BP)then
                       cycle
                    endif
                 endif
                 write (lunout, "(a7, 7(1xi6), 2(f12.4), (1xi6), (f12.4))", ADVANCE='NO') &
                      'contact', icon, iunit1, iunit2, imp1, imp2, &
                      imp1un, imp2un, &
                      go_nat(icon), factor_go(icon), &
                      icon_dummy_mgo(icon), coef_go(icon)
                 write(lunout, '(a4)') gotype2str()
              end if
           end do
           write (lunout, '(a)') ''
        end do
     end do
   
     write (lunout, '(a4)') '>>>>'
     write (lunout, '(a)') ''
  endif  ! ncon > 0


  ! ------------------------------------------------------------------
  ! LJ
  if (nLJ > 0) then
     write (lunout, '(a)') '<<<< LJ'
     write (lunout, '(a, i6)') '** total_contact = ', nLJ
     
     do iunit = 1, nunit_all
        do junit = iunit, nunit_all
           write (lunout, '(a, i6, a, i6)') '** contact between unit ', iunit, ' and ', junit
           !write (lunout, '(a, i6)') '** total_contact_unit = ', ncon_unit(iunit, junit) - nrna_bp_unit(iunit,junit)
           write (lunout, '(a, i6)') '** total_contact_unit = ', ncon_unit(iunit, junit) 
           write (lunout, '(a)', ADVANCE='NO') '**         iLJ iunit1-iunit2   imp1 - imp2 imp1un-imp2un'
           write (lunout, '(a)') '    distance        coef'
           do icon = 1, nLJ
              imp1 = iLJ2mp(1, icon)
              imp2 = iLJ2mp(2, icon)
              iunit1 = imp2unit(imp1)
              iunit2 = imp2unit(imp2)
              imp1un = imp1 - lunit2mp(1, iunit1) + 1
              imp2un = imp2 - lunit2mp(1, iunit2) + 1
                 
              if(iunit == iunit1 .and. junit == iunit2) then
                 if (iclass_unit(iunit) == CLASS%RNA .OR. iclass_unit(junit) == CLASS%RNA) then
                    if (icon2type(icon) == CONTYPE%RNA_BP)then
                       cycle
                    endif
                 endif
                 write (lunout, "(a7, 7(1xi6), 2(f12.4))") &
                      'contact', icon, iunit1, iunit2, imp1, imp2, &
                      imp1un, imp2un, &
                      LJ_nat(icon), coef_LJ(icon)
              end if
           end do
           !write (lunout, '(a)') ''
        end do
     end do
   
     write (lunout, '(a4)') '>>>>'
     write (lunout, '(a)') ''
  endif  ! nLJ > 0


!  ! ------------------------------------------------------------------
!  ! write the base-pair(RNA) interaction
!  if (nrna_bp > 0) then
!     write (lunout, '(a)') '<<<< native basepair '
!     write (lunout, '(a, i6)') '** total_contact = ', nrna_bp
!     write (lunout, '(a, f10.2, a)') '** definition_of_contact = ', inrna%dfcontact_bp, ' A'
!     write (lunout, '(a)') '** coef_go(kcal/mol) = factor_go * icon_dummy_mgo * cbp1210 * energy_unit_protein'
!     write (lunout, '(a)') ''
!     
!     do iunit = 1, nunit_all
!        do junit = iunit, nunit_all
!           if (iclass_unit(iunit) /= CLASS%RNA .OR. iclass_unit(junit) /= CLASS%RNA) then
!              cycle
!           endif
!           write (lunout, '(a, i6, a, i6)') '** contact between unit ', iunit, ' and ', junit
!           !write (lunout, '(a, i6)') '** total_contact_unit = ', ncon_unit(iunit, junit)
!           write (lunout, '(a)', ADVANCE='NO') '**        icon iunit1-iunit2   imp1 - imp2 imp1un-imp2un'
!           write (lunout, '(a)') '      go_nat   factor_go  dummy     coef_go'
!           do ibp = 1, nrna_bp
!              imp1 = irna_bp2mp(1, ibp)
!              imp2 = irna_bp2mp(2, ibp)
!              iunit1 = imp2unit(imp1)
!              iunit2 = imp2unit(imp2)
!              if(iunit == iunit1 .and. junit == iunit2) then
!                 imp1un = imp1 - lunit2mp(1, iunit1) + 1
!                 imp2un = imp2 - lunit2mp(1, iunit2) + 1
!                 
!                 if (iclass_unit(iunit) /= CLASS%RNA .OR. iclass_unit(junit) /= CLASS%RNA) then
!                    error_message = 'Error: logical defect in write_native_info (RNA_ST, but class != RNA)'
!                    call util_error(ERROR%STOP_ALL, error_message)
!                 endif
!   
!                 write (lunout, "(a8, 7(1xi6), 2(f12.4), (1xi6), (f12.4))", ADVANCE='NO') &
!                       'basepair', ibp, iunit1, iunit2, imp1, imp2, &
!                       imp1un, imp2un, rna_bp_nat(ibp), factor_rna_bp(ibp), &
!                       icon_dummy_mgo(ibp), coef_rna_bp(ibp)
!                 write (lunout, '(1xa4)', ADVANCE='NO') bptype2str()
!                 write (lunout, '(1xi1)') nhb_bp(ibp)
!              end if
!           end do
!           write (lunout, '(a)') ''
!        end do
!     end do
!   
!     write (lunout, '(a4)') '>>>>'
!     write (lunout, '(a)') ''
!  endif  ! nrna_bp > 0
!
!
!  ! ------------------------------------------------------------------
!  ! write the base-stack(RNA) interaction
!  if (nrna_st > 0) then
!     write (lunout, '(a)') '<<<< native basestack '
!     write (lunout, '(a, i6)') '** total_contact = ', nrna_st
!     write (lunout, '(a, f10.2, a)') '** definition_of_contact = ', inrna%dfcontact_st, ' A'
!     write (lunout, '(a)') '** coef_go(kcal/mol) = factor_go * icon_dummy_mgo * cst1210 * energy_unit_protein'
!     write (lunout, '(a)') ''
!     
!     do icon = 1, nrna_st
!        imp1 = irna_st2mp(1, icon)
!        imp2 = irna_st2mp(2, icon)
!        iunit1 = imp2unit(imp1)
!        iunit2 = imp2unit(imp2)
!        imp1un = imp1 - lunit2mp(1, iunit1) + 1
!        imp2un = imp2 - lunit2mp(1, iunit2) + 1
!                 
!        if (iunit1 /= iunit2) then
!           error_message = 'Error: logical defect in write_native_info (iunit1 /= iunit2 in RNA_ST)'
!           call util_error(ERROR%STOP_ALL, error_message)
!        endif
!   
!        if (iclass_unit(iunit1) /= CLASS%RNA) then
!           error_message = 'Error: logical defect in write_native_info (class != RNA in RNA_ST)'
!           call util_error(ERROR%STOP_ALL, error_message)
!        endif
!   
!        write (lunout, "(a9, 7(1xi6), 2(f12.4), (1xi6), (f12.4))", ADVANCE='NO') &
!              'basestack', icon, iunit1, iunit2, imp1, imp2, &
!              imp1un, imp2un, &
!              rna_st_nat(icon), factor_rna_st(icon), &
!              icon_dummy_mgo(icon), coef_rna_st(icon)
!        write(lunout, '(a4)') bstype2str()
!     end do
!   
!     write (lunout, '(a4)') '>>>>'
!  endif  ! nrna_st > 0

  ! ------------------------------------------------------------------
  ! write the base-stack(DT-RNA) interaction
  if (ndtrna_st > 0) then
     write (lunout, '(a)') '<<<< basestack of DT-RNA'
     write (lunout, '(a, i6)') '** total_contact = ', ndtrna_st
     write (lunout, '(a)') ''

     idih = 0
     do icon = 1, ndtrna_st
        imp1 = idtrna_st2mp(1,icon)
        imp2 = idtrna_st2mp(2,icon)
        iunit1 = imp2unit(imp1)
        iunit2 = imp2unit(imp2)
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit2) + 1
                 
        if (iunit1 /= iunit2) then
           error_message = 'Error: logical defect in write_native_info (iunit1 /= iunit2 in RNA_ST)'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
   
        if (iclass_unit(iunit1) /= CLASS%RNA) then
           error_message = 'Error: logical defect in write_native_info (class != RNA in RNA_ST)'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
   
        write (lunout, "(a7, 7(1xi6), 3(1xf12.4))")&
              'bs-dist', icon, iunit1, iunit2, imp1, imp2, imp1un, imp2un, &
              coef_dtrna_st(0,icon,IREP), dtrna_st_nat(1,icon), coef_dtrna_st(1,icon,IREP)
        !write(lunout, '(a4)') bstype2str()

        imp1 = idtrna_st2mp(3,icon)
        imp2 = idtrna_st2mp(4,icon)
        imp3 = idtrna_st2mp(5,icon)
        imp4 = idtrna_st2mp(6,icon)
        iunit1 = imp2unit(imp1)
        iunit2 = imp2unit(imp2)
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit1) + 1
        imp3un = imp3 - lunit2mp(1, iunit1) + 1
        imp4un = imp4 - lunit2mp(1, iunit1) + 1
        idih = idih + 1
        write (lunout, "(a7, 12(1xi6), 2(1xf12.4),1xa4)")&
             'bs-dihd', icon, idih, iunit1, iunit2, &
             imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
             dtrna_st_nat(2,icon) * 180.0e0_PREC / F_PI, coef_dtrna_st(2,icon,IREP), 'PSPS'

        imp1 = imp2
        imp2 = imp3
        imp3 = imp4
        imp4 = idtrna_st2mp(7,icon)
        imp1un = imp2un
        imp2un = imp3un
        imp3un = imp4un
        imp4un = imp4 - lunit2mp(1, iunit1) + 1
        idih = idih + 1
        write (lunout, "(a7, 12(1xi6), 2(1xf12.4),1xa4)")&
             'bs-dihd', icon, idih, iunit1, iunit2, &
             imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
             dtrna_st_nat(3,icon) * 180.0e0_PREC / F_PI, coef_dtrna_st(3,icon,IREP), 'SPSP'
     end do
   
     write (lunout, '(a4)') '>>>>'
  endif  ! ndtrna_st > 0

  ! ------------------------------------------------------------------
  ! write hydrogen bond (DT-RNA) interaction
  if (ndtrna_hb > 0) then
     write (lunout, '(a)') '<<<< hydrogen bond of DT-RNA'
     write (lunout, '(a, i6)') '** total_contact = ', ndtrna_hb
     write (lunout, '(a)') ''

     iba = 0
     idih = 0
     do icon = 1, ndtrna_hb
        imp1 = idtrna_hb2mp(1,icon)
        imp2 = idtrna_hb2mp(2,icon)
        iunit1 = imp2unit(imp1)
        iunit2 = imp2unit(imp2)
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit2) + 1
   
        if (iclass_unit(iunit1) /= CLASS%RNA) then
           error_message = 'Error: logical defect in write_native_info (class != RNA in RNA_ST)'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
   
        write (lunout, "(a7, 7(1xi6), 3(1xf12.4))")&
              'hb-dist', icon, iunit1, iunit2, imp1, imp2, imp1un, imp2un, &
              coef_dtrna_hb(0,icon), dtrna_hb_nat(1,icon), coef_dtrna_hb(1,icon)
        !write(lunout, '(a4)') bstype2str()

        imp1 = idtrna_hb2mp(3,icon)
        imp2 = idtrna_hb2mp(1,icon)
        imp3 = idtrna_hb2mp(2,icon)
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit1) + 1
        imp3un = imp3 - lunit2mp(1, iunit1) + 1
        iba = iba + 1
        write (lunout, "(a7, 10(1xi6), 2(1xf12.4))")&
             'hb-angl', icon, iba, iunit1, iunit2, &
             imp1, imp2, imp3, imp1un, imp2un, imp3un, &
             dtrna_hb_nat(2,icon) * 180.0e0_PREC / F_PI, coef_dtrna_hb(2,icon)

        imp1 = idtrna_hb2mp(1,icon)
        imp2 = idtrna_hb2mp(2,icon)
        imp3 = idtrna_hb2mp(4,icon)
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit1) + 1
        imp3un = imp3 - lunit2mp(1, iunit1) + 1
        iba = iba + 1
        write (lunout, "(a7, 10(1xi6), 2(1xf12.4))")&
             'hb-angl', icon, iba, iunit1, iunit2, &
             imp1, imp2, imp3, imp1un, imp2un, imp3un, &
             dtrna_hb_nat(3,icon) * 180.0e0_PREC / F_PI, coef_dtrna_hb(3,icon)

        imp1 = idtrna_hb2mp(3,icon)
        imp2 = idtrna_hb2mp(1,icon)
        imp3 = idtrna_hb2mp(2,icon)
        imp4 = idtrna_hb2mp(4,icon)
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit1) + 1
        imp3un = imp3 - lunit2mp(1, iunit1) + 1
        imp4un = imp4 - lunit2mp(1, iunit1) + 1
        idih = idih + 1
        write (lunout, "(a7, 12(1xi6), 2(1xf12.4),1xa4)")&
             'hb-dihd', icon, idih, iunit1, iunit2, &
             imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
             dtrna_hb_nat(4,icon) * 180.0e0_PREC / F_PI, coef_dtrna_hb(4,icon)

        imp1 = idtrna_hb2mp(5,icon)
        imp2 = idtrna_hb2mp(3,icon)
        imp3 = idtrna_hb2mp(1,icon)
        imp4 = idtrna_hb2mp(2,icon)
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit1) + 1
        imp3un = imp3 - lunit2mp(1, iunit1) + 1
        imp4un = imp4 - lunit2mp(1, iunit1) + 1
        idih = idih + 1
        write (lunout, "(a7, 12(1xi6), 2(1xf12.4),1xa4)")&
             'hb-dihd', icon, idih, iunit1, iunit2, &
             imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
             dtrna_hb_nat(5,icon) * 180.0e0_PREC / F_PI, coef_dtrna_hb(5,icon)

        imp1 = idtrna_hb2mp(1,icon)
        imp2 = idtrna_hb2mp(2,icon)
        imp3 = idtrna_hb2mp(4,icon)
        imp4 = idtrna_hb2mp(6,icon)
        imp1un = imp1 - lunit2mp(1, iunit1) + 1
        imp2un = imp2 - lunit2mp(1, iunit1) + 1
        imp3un = imp3 - lunit2mp(1, iunit1) + 1
        imp4un = imp4 - lunit2mp(1, iunit1) + 1
        idih = idih + 1
        write (lunout, "(a7, 12(1xi6), 2(1xf12.4),1xa4)")&
             'hb-dihd', icon, idih, iunit1, iunit2, &
             imp1, imp2, imp3, imp4, imp1un, imp2un, imp3un, imp4un, &
             dtrna_hb_nat(6,icon) * 180.0e0_PREC / F_PI, coef_dtrna_hb(6,icon)
     end do
   
     write (lunout, '(a4)') '>>>>'
  endif  ! ndtrna_st > 0
#ifdef MPI_PAR
  endif
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
   
   character(3) function bondtype2str()
      if (ibd2type(ibd) == BDTYPE%RNA_PS) then
         bondtype2str = ' PS'
      else if (ibd2type(ibd) == BDTYPE%RNA_SR .OR. ibd2type(ibd) == BDTYPE%RNA_SY) then
         if      (cmp2atom(imp2) == ' Ab ') then
            bondtype2str = ' SA'
         else if (cmp2atom(imp2) == ' Ub ') then
            bondtype2str = ' SU'
         else if (cmp2atom(imp2) == ' Gb ') then
            bondtype2str = ' SG'
         else if (cmp2atom(imp2) == ' Cb ') then
            bondtype2str = ' SC'
         else if (cmp2atom(imp2) == ' Rb ') then
            bondtype2str = ' SR'
         else if (cmp2atom(imp2) == ' Yb ') then
            bondtype2str = ' SY'
         else if (cmp2atom(imp2) == ' Nb ') then
            bondtype2str = ' SN'
         else
            error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) ' &
                           &//cmp2atom(imp2)
            write(*,*) imp2
            call util_error(ERROR%STOP_ALL, error_message)
         endif
      else if (ibd2type(ibd) == BDTYPE%RNA_SP) then
         bondtype2str = ' SP'
      else
         error_message = 'Error: logical defect in write_native_info (undefined BDTYPE)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   endfunction bondtype2str

   character(4) function angletype2str()
        if (iba2type(iba) == BATYPE%RNA_RSP .OR. iba2type(iba) == BATYPE%RNA_YSP) then
           if (cmp2atom(imp1) == ' Ab ') then
              angletype2str = ' ASP'
           else if (cmp2atom(imp1) == ' Ub ') then
              angletype2str = ' USP'
           else if (cmp2atom(imp1) == ' Gb ') then
              angletype2str = ' GSP'
           else if (cmp2atom(imp1) == ' Cb ') then
              angletype2str = ' CSP'
           else if (cmp2atom(imp1) == ' Rb ') then
              angletype2str = ' RSP'
           else if (cmp2atom(imp1) == ' Yb ') then
              angletype2str = ' YSP'
           else if (cmp2atom(imp1) == ' Nb ') then
              angletype2str = ' NSP'
           else
              error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                             &//cmp2atom(imp1)
               write(*,*) imp1
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        else if (iba2type(iba) == BATYPE%RNA_PSP) then
           angletype2str = ' PSP'
        else if (iba2type(iba) == BATYPE%RNA_SPS) then
           angletype2str = ' SPS'
        else if (iba2type(iba) == BATYPE%RNA_PSR .OR. iba2type(iba) == BATYPE%RNA_PSY) then
           if (cmp2atom(imp3) == ' Ab ') then
              angletype2str = ' PSA'
           else if (cmp2atom(imp3) == ' Ub ') then
              angletype2str = ' PSU'
           else if (cmp2atom(imp3) == ' Gb ') then
              angletype2str = ' PSG'
           else if (cmp2atom(imp3) == ' Cb ') then
              angletype2str = ' PSC'
           else if (cmp2atom(imp3) == ' Rb ') then
              angletype2str = ' PSR'
           else if (cmp2atom(imp3) == ' Yb ') then
              angletype2str = ' PSY'
           else if (cmp2atom(imp3) == ' Nb ') then
              angletype2str = ' PSN'
           else
              error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                             &//cmp2atom(imp3)
               write(*,*) imp3
              call util_error(ERROR%STOP_ALL, error_message)
           endif
!        else if (iba2type(iba) == BATYPE%RNA_BSB) then
!           angletype2str = ' BSB'
        else
           error_message = 'Error: logical defect in write_native_info (undefined BATYPE)'
           call util_error(ERROR%STOP_ALL, error_message)
        endif
   endfunction angletype2str

   character(5) function dihtype2str()
      if (idih2type(idih) == DIHTYPE%RNA_RSPS .OR. idih2type(idih) == DIHTYPE%RNA_YSPS) then
         if      (cmp2atom(imp1) == ' Ab ') then
            dihtype2str = ' ASPS'
         else if (cmp2atom(imp1) == ' Ub ') then
            dihtype2str = ' USPS'
         else if (cmp2atom(imp1) == ' Gb ') then
            dihtype2str = ' GSPS'
         else if (cmp2atom(imp1) == ' Cb ') then
            dihtype2str = ' CSPS'
         else if (cmp2atom(imp1) == ' Rb ') then
            dihtype2str = ' RSPS'
         else if (cmp2atom(imp1) == ' Yb ') then
            dihtype2str = ' YSPS'
         else if (cmp2atom(imp1) == ' Nb ') then
            dihtype2str = ' NSPS'
         else
            error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                           &//cmp2atom(imp1)
            call util_error(ERROR%STOP_ALL, error_message)
         endif
      else if (idih2type(idih) == DIHTYPE%RNA_PSPS) then
         dihtype2str = ' PSPS'
      else if (idih2type(idih) == DIHTYPE%RNA_SPSR .OR. idih2type(idih) == DIHTYPE%RNA_SPSY) then
         if      (cmp2atom(imp4) == ' Ab ') then
            dihtype2str = ' SPSA'
         else if (cmp2atom(imp4) == ' Ub ') then
            dihtype2str = ' SPSU'
         else if (cmp2atom(imp4) == ' Gb ') then
            dihtype2str = ' SPSG'
         else if (cmp2atom(imp4) == ' Cb ') then
            dihtype2str = ' SPSC'
         else if (cmp2atom(imp4) == ' Rb ') then
            dihtype2str = ' SPSR'
         else if (cmp2atom(imp4) == ' Yb ') then
            dihtype2str = ' SPSY'
         else if (cmp2atom(imp4) == ' Nb ') then
            dihtype2str = ' SPSN'
         else
            error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                           &//cmp2atom(imp4)
            call util_error(ERROR%STOP_ALL, error_message)
         endif
      else if (idih2type(idih) == DIHTYPE%RNA_SPSP) then
         dihtype2str = ' SPSP'
      else
         error_message = 'Error: logical defect in write_native_info (undefined DIHTYPE)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   endfunction dihtype2str
   
   character(4) function gotype2str()
      if (iclass_unit(iunit) == CLASS%PRO .AND. iclass_unit(junit) == CLASS%PRO) then
         gotype2str = ' p-p'
      else if (iclass_unit(iunit) == CLASS%RNA .OR. iclass_unit(junit) == CLASS%RNA) then
         if (icon2type(icon) == CONTYPE%PRO_RP) then
            gotype2str = ' p-P'
         else if (icon2type(icon) == CONTYPE%PRO_RS) then
            gotype2str = ' p-S'
         else if (icon2type(icon) == CONTYPE%PRO_RB) then
            gotype2str(1:3) = ' p-'
            if      (cmp2atom(imp1) == ' Ab ' .OR. cmp2atom(imp2) == ' Ab ') then
               gotype2str(4:4) = 'A'
            else if (cmp2atom(imp1) == ' Ub ' .OR. cmp2atom(imp2) == ' Ub ') then
               gotype2str(4:4) = 'U'
            else if (cmp2atom(imp1) == ' Gb ' .OR. cmp2atom(imp2) == ' Gb ') then
               gotype2str(4:4) = 'G'
            else if (cmp2atom(imp1) == ' Cb ' .OR. cmp2atom(imp2) == ' Cb ') then
               gotype2str(4:4) = 'C'
            else if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
               gotype2str(4:4) = 'R'
            else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
               gotype2str(4:4) = 'Y'
            else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
               gotype2str(4:4) = 'N'
            else
               error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                 &//cmp2atom(imp1)//','//cmp2atom(imp2)
               write(*,*) imp1,imp2
               call util_error(ERROR%STOP_ALL, error_message)
            endif
         else if (icon2type(icon) == CONTYPE%RP_RP) then
            gotype2str = ' P-P'
         else if (icon2type(icon) == CONTYPE%RP_RS) then
            gotype2str = ' P-S'
         else if (icon2type(icon) == CONTYPE%RP_RB) then
            gotype2str(1:3) = ' P-'
            if      (cmp2atom(imp1) == ' Ab ' .OR. cmp2atom(imp2) == ' Ab ') then
               gotype2str(4:4) = 'A'
            else if (cmp2atom(imp1) == ' Ub ' .OR. cmp2atom(imp2) == ' Ub ') then
               gotype2str(4:4) = 'U'
            else if (cmp2atom(imp1) == ' Gb ' .OR. cmp2atom(imp2) == ' Gb ') then
               gotype2str(4:4) = 'G'
            else if (cmp2atom(imp1) == ' Cb ' .OR. cmp2atom(imp2) == ' Cb ') then
               gotype2str(4:4) = 'C'
            else if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
               gotype2str(4:4) = 'R'
            else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
               gotype2str(4:4) = 'Y'
            else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
               gotype2str(4:4) = 'N'
            else
               error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                  &//cmp2atom(imp1)//','//cmp2atom(imp2)
               write(*,*) imp1,imp2
               call util_error(ERROR%STOP_ALL, error_message)
            endif
         else if (icon2type(icon) == CONTYPE%RS_RS) then
            gotype2str = ' S-S'
         else if (icon2type(icon) == CONTYPE%RS_RB) then
            gotype2str(1:3) = ' S-'
            if      (cmp2atom(imp1) == ' Ab ' .OR. cmp2atom(imp2) == ' Ab ') then
               gotype2str(4:4) = 'A'
            else if (cmp2atom(imp1) == ' Ub ' .OR. cmp2atom(imp2) == ' Ub ') then
               gotype2str(4:4) = 'U'
            else if (cmp2atom(imp1) == ' Gb ' .OR. cmp2atom(imp2) == ' Gb ') then
               gotype2str(4:4) = 'G'
            else if (cmp2atom(imp1) == ' Cb ' .OR. cmp2atom(imp2) == ' Cb ') then
               gotype2str(4:4) = 'C'
            else if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
               gotype2str(4:4) = 'R'
            else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
               gotype2str(4:4) = 'Y'
            else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
               gotype2str(4:4) = 'N'
            else
               error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                  &//cmp2atom(imp1)//','//cmp2atom(imp2)
               write(*,*) imp1,imp2
               call util_error(ERROR%STOP_ALL, error_message)
            endif
         else if (icon2type(icon) == CONTYPE%RB_RB) then
            if (cmp2atom(imp1) == ' Ab ' .OR. cmp2atom(imp2) == ' Ab ') then
               gotype2str(1:3) = ' A-'
               if      (cmp2atom(imp1) == ' Ub ' .OR. cmp2atom(imp2) == ' Ub ') then
                  gotype2str(4:4) = 'U'
               else if (cmp2atom(imp1) == ' Gb ' .OR. cmp2atom(imp2) == ' Gb ') then
                  gotype2str(4:4) = 'G'
               else if (cmp2atom(imp1) == ' Cb ' .OR. cmp2atom(imp2) == ' Cb ') then
                  gotype2str(4:4) = 'C'
               else if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
                  gotype2str(4:4) = 'R'
               else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
                  gotype2str(4:4) = 'Y'
               else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
                  gotype2str(4:4) = 'N'
               else if (cmp2atom(imp1) == ' Ab ' .AND. cmp2atom(imp2) == ' Ab ') then
                  gotype2str(4:4) = 'A'
               else
                  error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                     &//cmp2atom(imp1)//','//cmp2atom(imp2)
                  write(*,*) imp1,imp2
                  call util_error(ERROR%STOP_ALL, error_message)
               endif
            else if (cmp2atom(imp1) == ' Ub ' .OR. cmp2atom(imp2) == ' Ub ' ) then
               gotype2str(1:3) = ' U-'
               if      (cmp2atom(imp1) == ' Gb ' .OR. cmp2atom(imp2) == ' Gb ') then
                  gotype2str(4:4) = 'G'
               else if (cmp2atom(imp1) == ' Cb ' .OR. cmp2atom(imp2) == ' Cb ') then
                  gotype2str(4:4) = 'C'
               else if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
                  gotype2str(4:4) = 'R'
               else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
                  gotype2str(4:4) = 'Y'
               else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
                  gotype2str(4:4) = 'N'
               else if (cmp2atom(imp1) == ' Ub ' .AND. cmp2atom(imp2) == ' Ub ') then
                  gotype2str(4:4) = 'U'
               else
                  error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                     &//cmp2atom(imp1)//','//cmp2atom(imp2)
                  write(*,*) imp1,imp2
                  call util_error(ERROR%STOP_ALL, error_message)
               endif
            else if (cmp2atom(imp1) == ' Gb ' .OR. cmp2atom(imp2) == ' Gb ' ) then
               gotype2str(1:3) = ' G-'
               if      (cmp2atom(imp1) == ' Cb ' .OR. cmp2atom(imp2) == ' Cb ') then
                  gotype2str(4:4) = 'C'
               else if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
                  gotype2str(4:4) = 'R'
               else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
                  gotype2str(4:4) = 'Y'
               else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
                  gotype2str(4:4) = 'N'
               else if (cmp2atom(imp1) == ' Gb ' .AND. cmp2atom(imp2) == ' Gb ') then
                  gotype2str(4:4) = 'G'
               else
                  error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                     &//cmp2atom(imp1)//','//cmp2atom(imp2)
                  write(*,*) imp1,imp2
                  call util_error(ERROR%STOP_ALL, error_message)
               endif
            else if (cmp2atom(imp1) == ' Cb ' .OR. cmp2atom(imp2) == ' Cb ' ) then
               gotype2str(1:3) = ' C-'
               if      (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
                  gotype2str(4:4) = 'R'
               else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
                  gotype2str(4:4) = 'Y'
               else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
                  gotype2str(4:4) = 'N'
               else if (cmp2atom(imp1) == ' Cb ' .AND. cmp2atom(imp2) == ' Cb ') then
                  gotype2str(4:4) = 'C'
               else
                  error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                     &//cmp2atom(imp1)//','//cmp2atom(imp2)
                  write(*,*) imp1,imp2
                  call util_error(ERROR%STOP_ALL, error_message)
               endif
            else if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ' ) then
               gotype2str(1:3) = ' R-'
               if      (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
                  gotype2str(4:4) = 'Y'
               else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
                  gotype2str(4:4) = 'N'
               else if (cmp2atom(imp1) == ' Rb ' .AND. cmp2atom(imp2) == ' Rb ') then
                  gotype2str(4:4) = 'R'
               else
                  error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                     &//cmp2atom(imp1)//','//cmp2atom(imp2)
                  write(*,*) imp1,imp2
                  call util_error(ERROR%STOP_ALL, error_message)
               endif
            else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ' ) then
               gotype2str(1:3) = ' Y-'
               if      (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
                  gotype2str(4:4) = 'N'
               else if (cmp2atom(imp1) == ' Yb ' .AND. cmp2atom(imp2) == ' Yb ') then
                  gotype2str(4:4) = 'Y'
               else
                  error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                     &//cmp2atom(imp1)//','//cmp2atom(imp2)
                  write(*,*) imp1,imp2
                  call util_error(ERROR%STOP_ALL, error_message)
               endif
            else if (cmp2atom(imp1) == ' Nb ' .AND. cmp2atom(imp2) == ' Nb ' ) then
               gotype2str = ' N-N'
            else
               error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
                  &//','//cmp2atom(imp1)//','//cmp2atom(imp2)
               write(*,*) icon2type(icon)
               call util_error(ERROR%STOP_ALL, error_message)
            endif
         else
            error_message = 'Error: logical defect in write_native_info (undefined CONTYPE)'
            call util_error(ERROR%STOP_ALL, error_message)
         endif
      else
         gotype2str = '    '
      endif
   endfunction gotype2str
   
   character(4) function bptype2str()
      if (cmp2atom(imp1) == ' Ab ' .OR. cmp2atom(imp2) == ' Ab ') then
         bptype2str(1:3) = ' A-'
         if      (cmp2atom(imp1) == ' Ub ' .OR. cmp2atom(imp2) == ' Ub ') then
            bptype2str(4:4) = 'U'
         else if (cmp2atom(imp1) == ' Gb ' .OR. cmp2atom(imp2) == ' Gb ') then
            bptype2str(4:4) = 'G'
         else if (cmp2atom(imp1) == ' Cb ' .OR. cmp2atom(imp2) == ' Cb ') then
            bptype2str(4:4) = 'C'
         else if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
            bptype2str(4:4) = 'R'
         else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
            bptype2str(4:4) = 'Y'
         else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
            bptype2str(4:4) = 'N'
         else if (cmp2atom(imp1) == ' Ab ' .AND. cmp2atom(imp2) == ' Ab ') then
            bptype2str(4:4) = 'A'
         else
            error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
               &//cmp2atom(imp1)//','//cmp2atom(imp2)
            call util_error(ERROR%STOP_ALL, error_message)
         endif
      else if (cmp2atom(imp1) == ' Ub ' .OR. cmp2atom(imp2) == ' Ub ') then
         bptype2str(1:3) = ' U-'
         if      (cmp2atom(imp1) == ' Gb ' .OR. cmp2atom(imp2) == ' Gb ') then
            bptype2str(4:4) = 'G'
         else if (cmp2atom(imp1) == ' Cb ' .OR. cmp2atom(imp2) == ' Cb ') then
            bptype2str(4:4) = 'C'
         else if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
            bptype2str(4:4) = 'R'
         else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
            bptype2str(4:4) = 'Y'
         else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
            bptype2str(4:4) = 'N'
         else if (cmp2atom(imp1) == ' Ub ' .AND. cmp2atom(imp2) == ' Ub ') then
            bptype2str(4:4) = 'U'
         else
            error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
               &//cmp2atom(imp1)//','//cmp2atom(imp2)
            call util_error(ERROR%STOP_ALL, error_message)
         endif
      else if (cmp2atom(imp1) == ' Gb ' .OR. cmp2atom(imp2) == ' Gb ') then
         bptype2str(1:3) = ' G-'
         if      (cmp2atom(imp1) == ' Cb ' .OR. cmp2atom(imp2) == ' Cb ') then
            bptype2str(4:4) = 'C'
         else if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
            bptype2str(4:4) = 'R'
         else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
            bptype2str(4:4) = 'Y'
         else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
            bptype2str(4:4) = 'N'
         else if (cmp2atom(imp1) == ' Gb ' .AND. cmp2atom(imp2) == ' Gb ') then
            bptype2str(4:4) = 'G'
         else
            error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
               &//cmp2atom(imp1)//','//cmp2atom(imp2)
            call util_error(ERROR%STOP_ALL, error_message)
         endif
      else if (cmp2atom(imp1) == ' Cb ' .OR. cmp2atom(imp2) == ' Cb ') then
         bptype2str(1:3) = ' C-'
         if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
            bptype2str(4:4) = 'R'
         else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
            bptype2str(4:4) = 'Y'
         else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
            bptype2str(4:4) = 'N'
         else if (cmp2atom(imp1) == ' Cb ' .AND. cmp2atom(imp2) == ' Cb ') then
            bptype2str(4:4) = 'C'
         else
            error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
               &//cmp2atom(imp1)//','//cmp2atom(imp2)
            call util_error(ERROR%STOP_ALL, error_message)
         endif
      else if (cmp2atom(imp1) == ' Rb ' .OR. cmp2atom(imp2) == ' Rb ') then
         bptype2str(1:3) = ' R-'
         if      (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
            bptype2str(4:4) = 'Y'
         else if (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
            bptype2str(4:4) = 'N'
         else if (cmp2atom(imp1) == ' Rb ' .AND. cmp2atom(imp2) == ' Rb ') then
            bptype2str(4:4) = 'R'
         else
            error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
               &//cmp2atom(imp1)//','//cmp2atom(imp2)
            call util_error(ERROR%STOP_ALL, error_message)
         endif
      else if (cmp2atom(imp1) == ' Yb ' .OR. cmp2atom(imp2) == ' Yb ') then
         bptype2str(1:3) = ' Y-'
         if      (cmp2atom(imp1) == ' Nb ' .OR. cmp2atom(imp2) == ' Nb ') then
            bptype2str(4:4) = 'N'
         else if (cmp2atom(imp1) == ' Yb ' .AND. cmp2atom(imp2) == ' Yb ') then
            bptype2str(4:4) = 'Y'
         else
            error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
               &//cmp2atom(imp1)//','//cmp2atom(imp2)
            call util_error(ERROR%STOP_ALL, error_message)
         endif
      else if (cmp2atom(imp1) == ' Nb ' .AND. cmp2atom(imp2) == ' Nb ') then
            bptype2str = ' N-N'
      else
         error_message = 'Error: logical defect in write_native_info (undefined cmp2atom) '&
            &//cmp2atom(imp1)//','//cmp2atom(imp2)
         call util_error(ERROR%STOP_ALL, error_message)
      endif
      
   endfunction bptype2str

   character(4) function bstype2str()
      if (cmp2atom(imp1) == ' Ab ') then
         bstype2str(1:3) = ' A-'
      else if (cmp2atom(imp1) == ' Ub ') then
         bstype2str(1:3) = ' U-'
      else if (cmp2atom(imp1) == ' Gb ') then
         bstype2str(1:3) = ' G-'
      else if (cmp2atom(imp1) == ' Cb ') then
         bstype2str(1:3) = ' C-'
      else if (cmp2atom(imp1) == ' Rb ') then
         bstype2str(1:3) = ' R-'
      else if (cmp2atom(imp1) == ' Yb ') then
         bstype2str(1:3) = ' Y-'
      else if (cmp2atom(imp1) == ' Nb ') then
         bstype2str(1:3) = ' N-'
      else
         error_message = 'Error: logical defect in write_native_info (undefined STACK)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
      if (cmp2atom(imp2) == ' Ab ') then
         bstype2str(4:4) = 'A'
      else if (cmp2atom(imp2) == ' Ub ') then
         bstype2str(4:4) = 'U'
      else if (cmp2atom(imp2) == ' Gb ') then
         bstype2str(4:4) = 'G'
      else if (cmp2atom(imp2) == ' Cb ') then
         bstype2str(4:4) = 'C'
      else if (cmp2atom(imp2) == ' Rb ') then
         bstype2str(4:4) = 'R'
      else if (cmp2atom(imp2) == ' Yb ') then
         bstype2str(4:4) = 'Y'
      else if (cmp2atom(imp2) == ' Nb ') then
         bstype2str(4:4) = 'N'
      else
         error_message = 'Error: logical defect in write_native_info (undefined STACK)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   endfunction bstype2str
end subroutine write_nativeinfo
