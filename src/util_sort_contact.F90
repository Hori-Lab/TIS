! util_sort_contact
!> @brief Sort native contact

! ************************************************************************
subroutine util_sort_contact(ipost2pre_con)

  use const_maxsize
  use var_struct, only : nunit_all, nmp_all, ncon, &
                         icon2mp, lmp2con, icon2unit, icon2type, &
                         go_nat, go_nat2, factor_go, coef_go, &
                         icon_dummy_mgo, ncon_unit
  implicit none

  ! -----------------------------------------------------------------------
  ! native contact 
  integer, intent(out) :: ipost2pre_con(nmp_all*MXMPCON)
  ! intent(inout) :: icon2mp, icon2unit, go_nat, factor_go
  ! intent(out) :: lmp2con, ncon_unit

  ! -----------------------------------------------------------------------
  ! local variables
  integer :: i, imp, jmp
  integer :: icon, jcon, jcon1, jcon2, iunit, junit
  integer :: id, itest, istart, iend

  ! -----------------------------------------------------------------------
  do icon = 1, ncon
     ipost2pre_con(icon) = icon
  end do

  ! -----------------------------------------------------------------------
  ! sort icon2mp
  do icon = 1, ncon - 1
     itest = 0
     do jcon = 1, ncon - icon
        jcon1 = ncon - jcon
        jcon2 = ncon - jcon + 1
        if(icon2mp(1, jcon2) < icon2mp(1, jcon1)) then
           call iswapvar(ipost2pre_con(jcon1), ipost2pre_con(jcon2))
           call iswapvar(icon2mp(1, jcon1), icon2mp(1, jcon2))
           call iswapvar(icon2mp(2, jcon1), icon2mp(2, jcon2))
           call iswapvar(icon2unit(1, jcon1), icon2unit(1, jcon2))
           call iswapvar(icon2unit(2, jcon1), icon2unit(2, jcon2))
           call iswapvar(icon2type(jcon1), icon2type(jcon2))
           call rswapvar(go_nat(jcon1), go_nat(jcon2))
           call rswapvar(go_nat2(jcon1), go_nat2(jcon2))
           call rswapvar(factor_go(jcon1), factor_go(jcon2))
           call rswapvar(coef_go(jcon1), coef_go(jcon2))
           call iswapvar(icon_dummy_mgo(jcon1), icon_dummy_mgo(jcon2))
           itest = 1
        end if
     end do
     if(itest == 0) exit
  end do
      
  ! -----------------------------------------------------------------------
  ! calc lmp2con(imp)
  do i = 1, nmp_all
     lmp2con(i) = 0
  end do

  imp = 1
  do icon = 1, ncon
     jmp = icon2mp(1, icon)
     id = jmp - imp

     if(id == 0) then
        lmp2con(imp) = icon
     else if(id == 1) then
        imp = imp + 1
        lmp2con(imp) = icon
     else
        do i = 1, id - 1
           imp = imp + 1
           lmp2con(imp) = icon - 1
        end do
        imp = imp + 1
        lmp2con(imp) = icon
     end if
  end do

  ! -----------------------------------------------------------------------
  do imp = 1, nmp_all
     if(imp == 1) then
        istart = 1
     else
        istart = lmp2con(imp - 1) + 1
     endif
     iend = lmp2con(imp)
     do icon = istart, iend - 1
        do jcon = 1, iend - icon
           jcon1 = iend - jcon
           jcon2 = iend - jcon + 1
           if(icon2mp(2, jcon2) < icon2mp(2, jcon1)) then
              call iswapvar(ipost2pre_con(jcon1), ipost2pre_con(jcon2))
              call iswapvar(icon2mp(1, jcon1), icon2mp(1, jcon2))
              call iswapvar(icon2mp(2, jcon1), icon2mp(2, jcon2))
              call iswapvar(icon2unit(1, jcon1), icon2unit(1, jcon2))
              call iswapvar(icon2unit(2, jcon1), icon2unit(2, jcon2))
              call iswapvar(icon2type(jcon1), icon2type(jcon2))
              call rswapvar(go_nat(jcon1), go_nat(jcon2))
              call rswapvar(go_nat2(jcon1), go_nat2(jcon2))
              call rswapvar(factor_go(jcon1), factor_go(jcon2))
              call rswapvar(coef_go(jcon1), coef_go(jcon2))
              call iswapvar(icon_dummy_mgo(jcon1), icon_dummy_mgo(jcon2))
           end if
        end do
     end do
  end do

  ! -----------------------------------------------------------------------
  ! calc ncon_unit
  ncon_unit(1:nunit_all,1:nunit_all) = 0

  do icon = 1, ncon
     iunit = icon2unit(1, icon)
     junit = icon2unit(2, icon)
     ncon_unit(iunit, junit) = ncon_unit(iunit, junit) + 1 
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

end subroutine util_sort_contact
