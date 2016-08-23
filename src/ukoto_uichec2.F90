! ukoto_uichek2
!> @brief This subroutine is to read the string of characters &
!>        from the files (input, paramter, pdb).

!subroutine ukoto_uichec2(lunout, ctmp01, kcheck)
subroutine ukoto_uichec2(ctmp01, kcheck)

  use const_maxsize
  implicit none

  ! ----------------------------------------------------------------------
!  integer, intent(in) :: lunout
  character(CARRAY_MXCOLM), intent(in) :: ctmp01
  character(4), intent(inout) :: kcheck

  ! ----------------------------------------------------------------------
  ! local variables

  ! ----------------------------------------------------------------------
  if(ctmp01(1:1) == '*') then
     kcheck = 'COMM'
  else if(ctmp01(1:4) == '>>>>') then
     kcheck = 'END '
  else if(ctmp01(1:4) == '<<<<') then
     kcheck = 'ERR '
  else
     kcheck = 'DATA'
  end if

end subroutine ukoto_uichec2
