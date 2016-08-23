! util_sort_contact
!> @brief Sort native contact

! ************************************************************************
subroutine util_sort_rna_st(ipost2pre_rna_st)

  use const_maxsize
  use var_struct, only : nunit_all, nmp_all, nrna_st, &
                         irna_st2mp, lmp2rna_st, irna_st2unit, &
                         rna_st_nat, rna_st_nat2, factor_rna_st, &
                         coef_rna_st, coef_rna_st_a, coef_rna_st_fD, &
                         irna_st_dummy_mgo, &
                         nrna_st_unit
  implicit none

  ! -----------------------------------------------------------------------
  ! native contact 
  integer, intent(out), optional :: ipost2pre_rna_st(nmp_all*MXMPRNABP)

  ! -----------------------------------------------------------------------
  ! local variables
  integer :: i, imp, jmp
  integer :: icon, jcon, jcon1, jcon2, iunit, junit
  integer :: id, itest, istart, iend

  ! -----------------------------------------------------------------------
  if (present(ipost2pre_rna_st)) then
     do icon = 1, nrna_st
        ipost2pre_rna_st(icon) = icon
     end do
  endif

  ! -----------------------------------------------------------------------
  ! sort irna_st2mp
  do icon = 1, nrna_st - 1
     itest = 0
     do jcon = 1, nrna_st - icon
        jcon1 = nrna_st - jcon
        jcon2 = nrna_st - jcon + 1
        if(irna_st2mp(1, jcon2) < irna_st2mp(1, jcon1)) then
           if (present(ipost2pre_rna_st)) then
              call iswapvar(ipost2pre_rna_st(jcon1), ipost2pre_rna_st(jcon2))
           endif
           call iswapvar(irna_st2mp(1, jcon1), irna_st2mp(1, jcon2))
           call iswapvar(irna_st2mp(2, jcon1), irna_st2mp(2, jcon2))
           call iswapvar(irna_st2unit(1, jcon1), irna_st2unit(1, jcon2))
           call iswapvar(irna_st2unit(2, jcon1), irna_st2unit(2, jcon2))
           call rswapvar(rna_st_nat(jcon1), rna_st_nat(jcon2))
           call rswapvar(rna_st_nat2(jcon1), rna_st_nat2(jcon2))
           call rswapvar(factor_rna_st(jcon1), factor_rna_st(jcon2))
           call rswapvar(coef_rna_st(jcon1), coef_rna_st(jcon2))
           call rswapvar(coef_rna_st_a(jcon1), coef_rna_st_a(jcon2))
           call rswapvar(coef_rna_st_fD(jcon1), coef_rna_st_fD(jcon2))
           call iswapvar(irna_st_dummy_mgo(jcon1), irna_st_dummy_mgo(jcon2))
           itest = 1
        end if
     end do
     if(itest == 0) exit
  end do
      
  ! -----------------------------------------------------------------------
  ! calc lmp2con(imp)
  do i = 1, nmp_all
     lmp2rna_st(i) = 0
  end do

  imp = 1
  do icon = 1, nrna_st
     jmp = irna_st2mp(1, icon)
     id = jmp - imp

     if(id == 0) then
        lmp2rna_st(imp) = icon
     else if(id == 1) then
        imp = imp + 1
        lmp2rna_st(imp) = icon
     else
        do i = 1, id - 1
           imp = imp + 1
           lmp2rna_st(imp) = icon - 1
        end do
        imp = imp + 1
        lmp2rna_st(imp) = icon
     end if
  end do

  ! -----------------------------------------------------------------------
  do imp = 1, nmp_all
     if(imp == 1) then
        istart = 1
     else
        istart = lmp2rna_st(imp - 1) + 1
     endif
     iend = lmp2rna_st(imp)
     do icon = istart, iend - 1
        do jcon = 1, iend - icon
           jcon1 = iend - jcon
           jcon2 = iend - jcon + 1
           if(irna_st2mp(2, jcon2) < irna_st2mp(2, jcon1)) then
              if (present(ipost2pre_rna_st)) then
                 call iswapvar(ipost2pre_rna_st(jcon1), ipost2pre_rna_st(jcon2))
              endif
              call iswapvar(irna_st2mp(1, jcon1), irna_st2mp(1, jcon2))
              call iswapvar(irna_st2mp(2, jcon1), irna_st2mp(2, jcon2))
              call iswapvar(irna_st2unit(1, jcon1), irna_st2unit(1, jcon2))
              call iswapvar(irna_st2unit(2, jcon1), irna_st2unit(2, jcon2))
              call rswapvar(rna_st_nat(jcon1), rna_st_nat(jcon2))
              call rswapvar(rna_st_nat2(jcon1), rna_st_nat2(jcon2))
              call rswapvar(factor_rna_st(jcon1), factor_rna_st(jcon2))
              call rswapvar(coef_rna_st(jcon1), coef_rna_st(jcon2))
              call iswapvar(irna_st_dummy_mgo(jcon1), irna_st_dummy_mgo(jcon2))
           end if
        end do
     end do
  end do

  ! -----------------------------------------------------------------------
  ! calc nrna_st_unit
  do iunit = 1, nunit_all
     do junit = 1, nunit_all
        nrna_st_unit(iunit, junit) = 0
     end do
  end do

  do icon = 1, nrna_st
     iunit = irna_st2unit(1, icon)
     junit = irna_st2unit(2, icon)
     nrna_st_unit(iunit, junit) = nrna_st_unit(iunit, junit) + 1 
  end do

  ! *********************************************************************
contains

  ! *********************************************************************
  subroutine iswapvar(i1, i2)

    implicit none

    ! ------------------------------------------------------------------
    integer, intent(inout) :: i1, i2

    ! ------------------------------------------------------------------
    ! local variables
    integer :: itmp

    ! ------------------------------------------------------------------
    itmp = i1
    i1 = i2
    i2 = itmp

  end subroutine iswapvar

  ! *********************************************************************
  subroutine rswapvar(x1, x2)

    implicit none

    ! ------------------------------------------------------------------
    real(PREC), intent(inout) :: x1, x2

    ! ------------------------------------------------------------------
    ! local variables
    real(PREC) :: xtmp

    ! ------------------------------------------------------------------
    xtmp = x1
    x1 = x2
    x2 = xtmp

  end subroutine rswapvar

end subroutine util_sort_rna_st
