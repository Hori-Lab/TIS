!ifunc_char2int
!> @brief Converts one-letter alphabet to integer &
!>        number associated to the character.

!********************************************************************
integer function ifunc_char2int(char1)

  use const_maxsize
  use const_index

  implicit none

  ! ----------------------------------------------------------------
  character, intent(in) :: char1

  ! ----------------------------------------------------------------
  ! local variables
  character(CARRAY_MSG_ERROR) :: error_message
  
  ! ----------------------------------------------------------------

  ifunc_char2int  = 1
  if(char1 == 'a')then
     ifunc_char2int = 1
  else if(char1 == 'b')then
     ifunc_char2int = 2
  else if(char1 == 'c')then
     ifunc_char2int = 3
  else if(char1 == 'd')then
     ifunc_char2int = 4
  else if(char1 == 'e')then
     ifunc_char2int = 5
  else if(char1 == 'f')then
     ifunc_char2int = 6
  else if(char1 == 'g')then
     ifunc_char2int = 7
  else if(char1 == 'h')then
     ifunc_char2int = 8
  else if(char1 == 'i')then
     ifunc_char2int = 9
  else if(char1 == 'j')then
     ifunc_char2int = 10
  else if(char1 == 'k')then
     ifunc_char2int = 11
  else if(char1 == 'l')then
     ifunc_char2int = 12
  else if(char1 == 'm')then
     ifunc_char2int = 13
  else if(char1 == 'n')then
     ifunc_char2int = 14
  else if(char1 == 'o')then
     ifunc_char2int = 15
  else if(char1 == 'p')then
     ifunc_char2int = 16
  else if(char1 == 'q')then
     ifunc_char2int = 17
  else if(char1 == 'r')then
     ifunc_char2int = 18
  else if(char1 == 's')then
     ifunc_char2int = 19
  else if(char1 == 't')then
     ifunc_char2int = 20
  else if(char1 == 'u')then
     ifunc_char2int = 21
  else if(char1 == 'v')then
     ifunc_char2int = 22
  else if(char1 == 'w')then
     ifunc_char2int = 23
  else if(char1 == 'x')then
     ifunc_char2int = 24
  else if(char1 == 'y')then
     ifunc_char2int = 25
  else if(char1 == 'z')then
     ifunc_char2int = 26
  else
     error_message = 'Error: in ifunc_char2int there is no symbol such ' // char1
     call util_error(ERROR%STOP_ALL, error_message)
  end if

end function ifunc_char2int
