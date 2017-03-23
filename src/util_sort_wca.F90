! util_sort_contact
!> @brief Sort native contact

! ************************************************************************
subroutine util_sort_wca()

  use const_maxsize
  use var_struct, only : nmp_all, nwca, &
                         iwca2mp, lmp2wca, iwca2unit, &
                         wca_nat, wca_nat2, coef_wca
  implicit none

  ! -----------------------------------------------------------------------
  ! native contact 
  ! intent(inout) :: iwca2mp, iwca2unit, wca_nat, factor_go
  ! intent(out) :: lmp2wca, nwca_unit

  ! -----------------------------------------------------------------------
  ! local variables
  integer :: i, imp, jmp
  integer :: icon, jcon, jcon1, jcon2, iunit, junit
  integer :: id, itest, istart, iend

  ! -----------------------------------------------------------------------
  ! sort iwca2mp
  do icon = 1, nwca - 1
     itest = 0
     do jcon = 1, nwca - icon
        jcon1 = nwca - jcon
        jcon2 = nwca - jcon + 1
        if(iwca2mp(1, jcon2) < iwca2mp(1, jcon1)) then
           call iswapvar(iwca2mp(1, jcon1), iwca2mp(1, jcon2))
           call iswapvar(iwca2mp(2, jcon1), iwca2mp(2, jcon2))
           call iswapvar(iwca2unit(1, jcon1), iwca2unit(1, jcon2))
           call iswapvar(iwca2unit(2, jcon1), iwca2unit(2, jcon2))
           call rswapvar(wca_nat(jcon1), wca_nat(jcon2))
           call rswapvar(wca_nat2(jcon1), wca_nat2(jcon2))
           call rswapvar(coef_wca(jcon1,1), coef_wca(jcon2,1))
           call rswapvar(coef_wca(jcon1,2), coef_wca(jcon2,2))
           itest = 1
        end if
     end do
     if(itest == 0) exit
  end do
      
  ! -----------------------------------------------------------------------
  ! calc lmp2wca(imp)
  do i = 1, nmp_all
     lmp2wca(i) = 0
  end do

  imp = 1
  do icon = 1, nwca
     jmp = iwca2mp(1, icon)
     id = jmp - imp

     if(id == 0) then
        lmp2wca(imp) = icon
     else if(id == 1) then
        imp = imp + 1
        lmp2wca(imp) = icon
     else
        do i = 1, id - 1
           imp = imp + 1
           lmp2wca(imp) = icon - 1
        end do
        imp = imp + 1
        lmp2wca(imp) = icon
     end if
  end do

  ! -----------------------------------------------------------------------
  do imp = 1, nmp_all
     if(imp == 1) then
        istart = 1
     else
        istart = lmp2wca(imp - 1) + 1
     endif
     iend = lmp2wca(imp)
     do icon = istart, iend - 1
        do jcon = 1, iend - icon
           jcon1 = iend - jcon
           jcon2 = iend - jcon + 1
           if(iwca2mp(2, jcon2) < iwca2mp(2, jcon1)) then
              call iswapvar(iwca2mp(1, jcon1), iwca2mp(1, jcon2))
              call iswapvar(iwca2mp(2, jcon1), iwca2mp(2, jcon2))
              call iswapvar(iwca2unit(1, jcon1), iwca2unit(1, jcon2))
              call iswapvar(iwca2unit(2, jcon1), iwca2unit(2, jcon2))
              call rswapvar(wca_nat(jcon1), wca_nat(jcon2))
              call rswapvar(wca_nat2(jcon1), wca_nat2(jcon2))
              call rswapvar(coef_wca(jcon1,1), coef_wca(jcon2,1))
              call rswapvar(coef_wca(jcon1,2), coef_wca(jcon2,2))
           end if
        end do
     end do
  end do


  do icon = 1, nwca
     iunit = iwca2unit(1, icon)
     junit = iwca2unit(2, icon)
     !nwca_unit(iunit, junit) = nwca_unit(iunit, junit) + 1 
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

end subroutine util_sort_wca
