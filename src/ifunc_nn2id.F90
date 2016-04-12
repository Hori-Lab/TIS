! ifunc_nn2id
!> @brief This function converts strings of nearest-neighbor nucleic acids &
!> into the identification number (id).
!********************************************************************
integer function ifunc_nn2id(str)

  use const_maxsize
  use const_index

  implicit none
  ! ----------------------------------------------------------------
  character(2), intent(in) :: str

  ! ----------------------------------------------------------------
  ! local variables
  character(CARRAY_MSG_ERROR) :: error_message
  
  ifunc_nn2id = 0
  if(str == 'AA')then
     ifunc_nn2id = NN%AA
  else if(str == 'AU') then
     ifunc_nn2id = NN%AU
  else if(str == 'AG') then
     ifunc_nn2id = NN%AG
  else if(str == 'AC') then
     ifunc_nn2id = NN%AC
  else if(str == 'UA') then
     ifunc_nn2id = NN%UA
  else if(str == 'UU') then
     ifunc_nn2id = NN%UU
  else if(str == 'UG') then
     ifunc_nn2id = NN%UG
  else if(str == 'UC') then
     ifunc_nn2id = NN%UC
  else if(str == 'GA') then
     ifunc_nn2id = NN%GA
  else if(str == 'GU') then
     ifunc_nn2id = NN%GU
  else if(str == 'GG') then
     ifunc_nn2id = NN%GG
  else if(str == 'GC') then
     ifunc_nn2id = NN%GC
  else if(str == 'CA') then
     ifunc_nn2id = NN%CA
  else if(str == 'CU') then
     ifunc_nn2id = NN%CU
  else if(str == 'CG') then
     ifunc_nn2id = NN%CG
  else if(str == 'CC') then
     ifunc_nn2id = NN%CC
  else
     error_message = 'Error: in ifunc_nn2id there is no name such ' // str
     call util_error(ERROR%STOP_ALL, error_message)
  end if

end function ifunc_nn2id
