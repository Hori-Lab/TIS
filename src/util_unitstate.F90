! util_unitstate
!> @brief Read a formatted string, 'u1-u2s', and store the & 
!>        information of unit numbers and state. 'u1' and 'u2' &
!>        are integers, representing the start and the end of &
!>        unit number, and 's' is an alphabet indicating the state

subroutine util_unitstate( char12, & ! [i ]
                           inunit, & ! [ o]
                           instate & ! [ o]
                          )
  use const_maxsize
  use const_index

  implicit none

  ! --------------------------------------------------------------------
  character(12), intent(in)  :: char12
  integer,       intent(out) :: inunit(2), instate

  ! --------------------------------------------------------------------
  ! local variables
  integer :: i, ie, iend(2)
  integer :: ifunc_char2int
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------

  inunit(1) = 0
  inunit(2) = 0
  instate = 0

  ie = 1
  iend(1) = 0
  iend(2) = 0
  do i = 1, 12
     if(char12(i:i) == ' ') then
        exit
     else if(char12(i:i) == '-') then
        iend(1) = i - 1
        iend(2) = i
        ie = ie + 1
     else if(char12(i:i) >= '0' .and. char12(i:i) <= '9') then
        iend(ie) = iend(ie) + 1
     else if(char12(i:i) >= 'a' .and. char12(i:i) <= 'f') then
        instate = ifunc_char2int(char12(i:i))
        exit
     else
        error_message = 'Error: invalid argument in state and status in util_unitstate'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
  end do

  read (char12(1:iend(1)), *) inunit(1)
  if(ie == 1) then
     inunit(2) = inunit(1)
  else
     read (char12(iend(1)+2:iend(2)), *) inunit(2)
  end if

end subroutine util_unitstate
