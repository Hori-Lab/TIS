! util_char_type
!> @brief check character whether number or not

subroutine util_char_type(char1, flag_num)

  implicit none

  ! --------------------------------------------------------------------
  character, intent(in) :: char1
  logical, intent(out) :: flag_num

  ! --------------------------------------------------------------------

  if(char1 == '1') then
     flag_num = .TRUE.
  else if(char1 == '2') then
     flag_num = .TRUE.
  else if(char1 == '3') then
     flag_num = .TRUE.
  else if(char1 == '4') then
     flag_num = .TRUE.
  else if(char1 == '5') then
     flag_num = .TRUE.
  else if(char1 == '6') then
     flag_num = .TRUE.
  else if(char1 == '7') then
     flag_num = .TRUE.
  else if(char1 == '8') then
     flag_num = .TRUE.
  else if(char1 == '9') then
     flag_num = .TRUE.
  else if(char1 == '0') then
     flag_num = .TRUE.
  else
     flag_num = .FALSE.
  end if

end subroutine util_char_type
