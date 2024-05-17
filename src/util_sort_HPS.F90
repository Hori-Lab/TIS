! util_sort_contact
!> @brief Sort native contact

! ************************************************************************
subroutine util_sort_HPS()

  use const_maxsize
  use var_struct, only : nmp_all, nHPS, &
                         iHPS2mp, lmp2HPS, iHPS2unit, &
                         HPS_nat, HPS_nat2, coef_HPS, lambda
  implicit none

  ! -----------------------------------------------------------------------
  ! native contact 
  ! intent(inout) :: iHPS2mp, iHPS2unit, HPS_nat, factor_go
  ! intent(out) :: lmp2HPS, nHPS_unit

  ! -----------------------------------------------------------------------
  ! local variables
  integer :: i, imp, jmp
  integer :: icon, jcon, jcon1, jcon2, iunit, junit
  integer :: id, itest, istart, iend

  ! -----------------------------------------------------------------------
  ! sort iHPS2mp
  do icon = 1, nHPS - 1
     itest = 0
     do jcon = 1, nHPS - icon
        jcon1 = nHPS - jcon
        jcon2 = nHPS - jcon + 1
        if(iHPS2mp(1, jcon2) < iHPS2mp(1, jcon1)) then
           call iswapvar(iHPS2mp(1, jcon1), iHPS2mp(1, jcon2))
           call iswapvar(iHPS2mp(2, jcon1), iHPS2mp(2, jcon2))
           call iswapvar(iHPS2unit(1, jcon1), iHPS2unit(1, jcon2))
           call iswapvar(iHPS2unit(2, jcon1), iHPS2unit(2, jcon2))
           call rswapvar(HPS_nat(jcon1), HPS_nat(jcon2))
           call rswapvar(HPS_nat2(jcon1), HPS_nat2(jcon2))
           call rswapvar(coef_HPS(jcon1), coef_HPS(jcon2))
           call rswapvar(lambda(jcon1), lambda(jcon2))
           itest = 1
        end if
     end do
     if(itest == 0) exit
  end do
      
  ! -----------------------------------------------------------------------
  ! calc lmp2HPS(imp)
  do i = 1, nmp_all
     lmp2HPS(i) = 0
  end do

  imp = 1
  do icon = 1, nHPS
     jmp = iHPS2mp(1, icon)
     id = jmp - imp

     if(id == 0) then
        lmp2HPS(imp) = icon
     else if(id == 1) then
        imp = imp + 1
        lmp2HPS(imp) = icon
     else
        do i = 1, id - 1
           imp = imp + 1
           lmp2HPS(imp) = icon - 1
        end do
        imp = imp + 1
        lmp2HPS(imp) = icon
     end if
  end do

  ! -----------------------------------------------------------------------
  do imp = 1, nmp_all
     if(imp == 1) then
        istart = 1
     else
        istart = lmp2HPS(imp - 1) + 1
     endif
     iend = lmp2HPS(imp)
     do icon = istart, iend - 1
        do jcon = 1, iend - icon
           jcon1 = iend - jcon
           jcon2 = iend - jcon + 1
           if(iHPS2mp(2, jcon2) < iHPS2mp(2, jcon1)) then
              call iswapvar(iHPS2mp(1, jcon1), iHPS2mp(1, jcon2))
              call iswapvar(iHPS2mp(2, jcon1), iHPS2mp(2, jcon2))
              call iswapvar(iHPS2unit(1, jcon1), iHPS2unit(1, jcon2))
              call iswapvar(iHPS2unit(2, jcon1), iHPS2unit(2, jcon2))
              call rswapvar(HPS_nat(jcon1), HPS_nat(jcon2))
              call rswapvar(HPS_nat2(jcon1), HPS_nat2(jcon2))
              call rswapvar(coef_HPS(jcon1), coef_HPS(jcon2))
           end if
        end do
     end do
  end do


  do icon = 1, nHPS
     iunit = iHPS2unit(1, icon)
     junit = iHPS2unit(2, icon)
     !nHPS_unit(iunit, junit) = nHPS_unit(iunit, junit) + 1 
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

end subroutine util_sort_HPS
