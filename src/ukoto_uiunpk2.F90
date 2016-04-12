! ukoto_uiunpk2
!> @brief This subroutine is to read the string of characters &
!>        from the files (input, paramter, pdb).

subroutine ukoto_uiunpk2(lunout, ctmp01, ifirst)

  use const_maxsize
  use const_index
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lunout
  integer, intent(out) :: ifirst
  character(CARRAY_MXCOLM), intent(in) :: ctmp01

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
  do i = 5, CARRAY_MXCOLM
     ifirst = i
     if(ctmp01(i:i) /= ' ') then
        exit
     end if
  end do

  if(ifirst > CARRAY_MXCOLM - 15) then
     error_message = 'Error: input field name is wrong.'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

end subroutine ukoto_uiunpk2
