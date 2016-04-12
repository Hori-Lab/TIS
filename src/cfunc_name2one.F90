! cfunc_name2one
!> @brief This function converts the amino-acid name from 3-letter to 1-letter.


! ********************************************************************
! Mar.10.2000
! ********************************************************************
character(1) function cfunc_name2one(name)

  use const_maxsize
  use const_index
  implicit none

!  character(1) :: cfunc_name2one

  ! -----------------------------------------------------------------
  character(3) :: name
  character(CARRAY_MSG_ERROR) :: error_message

  ! -----------------------------------------------------------------

  if(name == 'ALA') then
     cfunc_name2one = 'A'
  else if(name == 'ARG') then
     cfunc_name2one = 'R'
  else if(name == 'ASN') then
     cfunc_name2one = 'N'
  else if(name == 'ASP') then
     cfunc_name2one = 'D'
  else if(name == 'CYS') then
     cfunc_name2one = 'C'
  else if(name == 'GLN') then
     cfunc_name2one = 'Q'
  else if(name == 'GLU') then
     cfunc_name2one = 'E'
  else if(name == 'GLY') then
     cfunc_name2one = 'G'
  else if(name == 'HIS' .or. name=='HIE' .or. name=='HID' .or. name=='HIP' .or. &
                             name=='HSE' .or. name=='HSD' .or. name=='HSP' ) then
     cfunc_name2one = 'H'
  else if(name == 'ILE') then
     cfunc_name2one = 'I'
  else if(name == 'LEU') then
     cfunc_name2one = 'L'
  else if(name == 'LYS') then
     cfunc_name2one = 'K'
  else if(name == 'MET') then
     cfunc_name2one = 'M'
  else if(name == 'PHE') then
     cfunc_name2one = 'F'
  else if(name == 'PRO') then
     cfunc_name2one = 'P'
  else if(name == 'SER') then
     cfunc_name2one = 'S'
  else if(name == 'THR') then
     cfunc_name2one = 'T'
  else if(name == 'TRP') then
     cfunc_name2one = 'W'
  else if(name == 'TYR') then
     cfunc_name2one = 'Y'
  else if(name == 'VAL') then
     cfunc_name2one = 'V'
  else if(name == ' DA') then
     cfunc_name2one = 'A'
  else if(name == ' DT') then
     cfunc_name2one = 'T'
  else if(name == ' DG') then
     cfunc_name2one = 'G'
  else if(name == ' DC') then
     cfunc_name2one = 'C'
  else if(name == 'COR') then
     cfunc_name2one = 'C'
  else if(name == 'INT') then
     cfunc_name2one = 'I'
  else if(name == 'TAI') then
     cfunc_name2one = 'T'
  else
     error_message = 'Error: in cfunc_name2one' // name
     call util_error(ERROR%STOP_ALL, error_message)
  end if

end function cfunc_name2one

