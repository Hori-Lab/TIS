! cfunc_int2char
!> @brief Map from integers (1-26) to alphabets (a-z)

character(1) function cfunc_int2char(i1)

  use const_maxsize
  use const_index

  implicit none

  ! ----------------------------------------------------------------
  integer, intent(in) :: i1

  ! ----------------------------------------------------------------
  ! local variables
  integer :: i1mod
  character(CARRAY_MSG_ERROR) :: error_message
  
  ! ----------------------------------------------------------------

  if (i1 == 0) then
     cfunc_int2char = ' '
     return
  endif

  i1mod = mod(i1 - 1, 26) + 1

  if(i1mod == 1) then
     cfunc_int2char = 'a'
  else if(i1mod == 2) then
     cfunc_int2char = 'b'
  else if(i1mod == 3) then
     cfunc_int2char = 'c'
  else if(i1mod == 4) then
     cfunc_int2char = 'd'
  else if(i1mod == 5) then
     cfunc_int2char = 'e'
  else if(i1mod == 6) then
     cfunc_int2char = 'f'
  else if(i1mod == 7) then
     cfunc_int2char = 'g'
  else if(i1mod == 8) then
     cfunc_int2char = 'h'
  else if(i1mod == 9) then
     cfunc_int2char = 'i'
  else if(i1mod == 10) then
     cfunc_int2char = 'j'
  else if(i1mod == 11) then
     cfunc_int2char = 'k'
  else if(i1mod == 12) then
     cfunc_int2char = 'l'
  else if(i1mod == 13) then
     cfunc_int2char = 'm'
  else if(i1mod == 14) then
     cfunc_int2char = 'n'
  else if(i1mod == 15) then
     cfunc_int2char = 'o'
  else if(i1mod == 16) then
     cfunc_int2char = 'p'
  else if(i1mod == 17) then
     cfunc_int2char = 'q'
  else if(i1mod == 18) then
     cfunc_int2char = 'r'
  else if(i1mod == 19) then
     cfunc_int2char = 's'
  else if(i1mod == 20) then
     cfunc_int2char = 't'
  else if(i1mod == 21) then
     cfunc_int2char = 'u'
  else if(i1mod == 22) then
     cfunc_int2char = 'v'
  else if(i1mod == 23) then
     cfunc_int2char = 'w'
  else if(i1mod == 24) then
     cfunc_int2char = 'x'
  else if(i1mod == 25) then
     cfunc_int2char = 'y'
  else if(i1mod == 26) then
     cfunc_int2char = 'z'
  else
     error_message = 'Error: in cfunc_int2char'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

end function cfunc_int2char
