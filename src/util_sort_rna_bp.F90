! util_sort_contact
!> @brief Sort native contact

! ************************************************************************
subroutine util_sort_rna_bp(ipost2pre_rna_bp)

  use const_maxsize
  use var_struct, only : nunit_all, nmp_all, nrna_bp, &
                         irna_bp2mp, lmp2rna_bp, irna_bp2unit, &
                         rna_bp_nat, rna_bp_nat2, factor_rna_bp, &
                         coef_rna_bp, coef_rna_bp_a, coef_rna_bp_fD, &
                         irna_bp_dummy_mgo, nhb_bp, &
                         nrna_bp_unit
  implicit none

  ! -----------------------------------------------------------------------
  ! native contact 
  integer, intent(out), optional :: ipost2pre_rna_bp(nmp_all*MXMPRNABP)

  ! -----------------------------------------------------------------------
  ! local variables
  integer :: i, imp, jmp
  integer :: icon, jcon, jcon1, jcon2, iunit, junit
  integer :: id, itest, istart, iend

  ! -----------------------------------------------------------------------
  if (present(ipost2pre_rna_bp)) then
     do icon = 1, nrna_bp
        ipost2pre_rna_bp(icon) = icon
     end do
  endif

  ! -----------------------------------------------------------------------
  ! sort irna_bp2mp
  do icon = 1, nrna_bp - 1
     itest = 0
     do jcon = 1, nrna_bp - icon
        jcon1 = nrna_bp - jcon
        jcon2 = nrna_bp - jcon + 1
        if(irna_bp2mp(1, jcon2) < irna_bp2mp(1, jcon1)) then
           if (present(ipost2pre_rna_bp)) then
              call iswapvar(ipost2pre_rna_bp(jcon1), ipost2pre_rna_bp(jcon2))
           endif
           call iswapvar(irna_bp2mp(1, jcon1), irna_bp2mp(1, jcon2))
           call iswapvar(irna_bp2mp(2, jcon1), irna_bp2mp(2, jcon2))
           call iswapvar(irna_bp2unit(1, jcon1), irna_bp2unit(1, jcon2))
           call iswapvar(irna_bp2unit(2, jcon1), irna_bp2unit(2, jcon2))
           call rswapvar(rna_bp_nat(jcon1), rna_bp_nat(jcon2))
           call rswapvar(rna_bp_nat2(jcon1), rna_bp_nat2(jcon2))
           call rswapvar(factor_rna_bp(jcon1), factor_rna_bp(jcon2))
           call rswapvar(coef_rna_bp(jcon1), coef_rna_bp(jcon2))
           call rswapvar(coef_rna_bp_a(jcon1), coef_rna_bp_a(jcon2))
           call rswapvar(coef_rna_bp_fD(jcon1), coef_rna_bp_fD(jcon2))
           call iswapvar(nhb_bp(jcon1), nhb_bp(jcon2))
           call iswapvar(irna_bp_dummy_mgo(jcon1), irna_bp_dummy_mgo(jcon2))
           itest = 1
        end if
     end do
     if(itest == 0) exit
  end do
      
  ! -----------------------------------------------------------------------
  ! calc lmp2con(imp)
  do i = 1, nmp_all
     lmp2rna_bp(i) = 0
  end do

  imp = 1
  do icon = 1, nrna_bp
     jmp = irna_bp2mp(1, icon)
     id = jmp - imp

     if(id == 0) then
        lmp2rna_bp(imp) = icon
     else if(id == 1) then
        imp = imp + 1
        lmp2rna_bp(imp) = icon
     else
        do i = 1, id - 1
           imp = imp + 1
           lmp2rna_bp(imp) = icon - 1
        end do
        imp = imp + 1
        lmp2rna_bp(imp) = icon
     end if
  end do

  ! -----------------------------------------------------------------------
  do imp = 1, nmp_all
     if(imp == 1) then
        istart = 1
     else
        istart = lmp2rna_bp(imp - 1) + 1
     endif
     iend = lmp2rna_bp(imp)
     do icon = istart, iend - 1
        do jcon = 1, iend - icon
           jcon1 = iend - jcon
           jcon2 = iend - jcon + 1
           if(irna_bp2mp(2, jcon2) < irna_bp2mp(2, jcon1)) then
              if (present(ipost2pre_rna_bp)) then
                 call iswapvar(ipost2pre_rna_bp(jcon1), ipost2pre_rna_bp(jcon2))
              endif
              call iswapvar(irna_bp2mp(1, jcon1), irna_bp2mp(1, jcon2))
              call iswapvar(irna_bp2mp(2, jcon1), irna_bp2mp(2, jcon2))
              call iswapvar(irna_bp2unit(1, jcon1), irna_bp2unit(1, jcon2))
              call iswapvar(irna_bp2unit(2, jcon1), irna_bp2unit(2, jcon2))
              call rswapvar(rna_bp_nat(jcon1), rna_bp_nat(jcon2))
              call rswapvar(rna_bp_nat2(jcon1), rna_bp_nat2(jcon2))
              call rswapvar(factor_rna_bp(jcon1), factor_rna_bp(jcon2))
              call rswapvar(coef_rna_bp(jcon1), coef_rna_bp(jcon2))
              call rswapvar(coef_rna_bp_a(jcon1), coef_rna_bp_a(jcon2))
              call rswapvar(coef_rna_bp_fD(jcon1), coef_rna_bp_fD(jcon2))
              call iswapvar(nhb_bp(jcon1), nhb_bp(jcon2))
              call iswapvar(irna_bp_dummy_mgo(jcon1), irna_bp_dummy_mgo(jcon2))
           end if
        end do
     end do
  end do

  ! -----------------------------------------------------------------------
  ! calc nrna_bp_unit
  do iunit = 1, nunit_all
     do junit = 1, nunit_all
        nrna_bp_unit(iunit, junit) = 0
     end do
  end do

  do icon = 1, nrna_bp
     iunit = irna_bp2unit(1, icon)
     junit = irna_bp2unit(2, icon)
     nrna_bp_unit(iunit, junit) = nrna_bp_unit(iunit, junit) + 1 
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

end subroutine util_sort_rna_bp
