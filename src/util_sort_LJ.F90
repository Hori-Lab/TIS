! util_sort_contact
!> @brief Sort native contact

! ************************************************************************
subroutine util_sort_LJ()

  use const_maxsize
  use var_struct, only : nmp_all, nLJ, &
                         iLJ2mp, lmp2LJ, iLJ2unit, &
                         LJ_nat, LJ_nat2, coef_LJ
  implicit none

  ! -----------------------------------------------------------------------
  ! native contact 
  ! intent(inout) :: iLJ2mp, iLJ2unit, LJ_nat, factor_go
  ! intent(out) :: lmp2LJ, nLJ_unit

  ! -----------------------------------------------------------------------
  ! local variables
  integer :: i, imp, jmp
  integer :: icon, jcon, jcon1, jcon2, iunit, junit
  integer :: id, itest, istart, iend

  ! -----------------------------------------------------------------------
  ! sort iLJ2mp
  do icon = 1, nLJ - 1
     itest = 0
     do jcon = 1, nLJ - icon
        jcon1 = nLJ - jcon
        jcon2 = nLJ - jcon + 1
        if(iLJ2mp(1, jcon2) < iLJ2mp(1, jcon1)) then
           call iswapvar(iLJ2mp(1, jcon1), iLJ2mp(1, jcon2))
           call iswapvar(iLJ2mp(2, jcon1), iLJ2mp(2, jcon2))
           call iswapvar(iLJ2unit(1, jcon1), iLJ2unit(1, jcon2))
           call iswapvar(iLJ2unit(2, jcon1), iLJ2unit(2, jcon2))
           call rswapvar(LJ_nat(jcon1), LJ_nat(jcon2))
           call rswapvar(LJ_nat2(jcon1), LJ_nat2(jcon2))
           call rswapvar(coef_LJ(jcon1), coef_LJ(jcon2))
           itest = 1
        end if
     end do
     if(itest == 0) exit
  end do
      
  ! -----------------------------------------------------------------------
  ! calc lmp2LJ(imp)
  do i = 1, nmp_all
     lmp2LJ(i) = 0
  end do

  imp = 1
  do icon = 1, nLJ
     jmp = iLJ2mp(1, icon)
     id = jmp - imp

     if(id == 0) then
        lmp2LJ(imp) = icon
     else if(id == 1) then
        imp = imp + 1
        lmp2LJ(imp) = icon
     else
        do i = 1, id - 1
           imp = imp + 1
           lmp2LJ(imp) = icon - 1
        end do
        imp = imp + 1
        lmp2LJ(imp) = icon
     end if
  end do

  ! -----------------------------------------------------------------------
  do imp = 1, nmp_all
     if(imp == 1) then
        istart = 1
     else
        istart = lmp2LJ(imp - 1) + 1
     endif
     iend = lmp2LJ(imp)
     do icon = istart, iend - 1
        do jcon = 1, iend - icon
           jcon1 = iend - jcon
           jcon2 = iend - jcon + 1
           if(iLJ2mp(2, jcon2) < iLJ2mp(2, jcon1)) then
              call iswapvar(iLJ2mp(1, jcon1), iLJ2mp(1, jcon2))
              call iswapvar(iLJ2mp(2, jcon1), iLJ2mp(2, jcon2))
              call iswapvar(iLJ2unit(1, jcon1), iLJ2unit(1, jcon2))
              call iswapvar(iLJ2unit(2, jcon1), iLJ2unit(2, jcon2))
              call rswapvar(LJ_nat(jcon1), LJ_nat(jcon2))
              call rswapvar(LJ_nat2(jcon1), LJ_nat2(jcon2))
              call rswapvar(coef_LJ(jcon1), coef_LJ(jcon2))
           end if
        end do
     end do
  end do


  do icon = 1, nLJ
     iunit = iLJ2unit(1, icon)
     junit = iLJ2unit(2, icon)
     !nLJ_unit(iunit, junit) = nLJ_unit(iunit, junit) + 1 
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

end subroutine util_sort_LJ
