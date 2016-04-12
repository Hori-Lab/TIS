! ukoto_uvvalue2
!> @brief This routine is for reading input data of 'koto' enclosed in 
!> '<<<< (naminp)' and '>>>>' excluding comments, where '(naminp)' can &
!> be any text
subroutine ukoto_rvalue2(lunout, cwork, rvalue, cvalue)

  use const_maxsize
  use const_index
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lunout
  real(PREC), intent(out) :: rvalue
  character(CARRAY_MXCOLM), intent(in) :: cwork(2)
  character(CARRAY_MXCOLM), intent(in) :: cvalue

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i !,len
  integer :: intrea
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
!  len = len_trim(cwork(1))

!  if(cwork(1)(1:len) /= cvalue(1:len)) then
  if(cwork(1) /= cvalue) then
     return
  end if

  !---- store characters into 'ctmp01' -----------------------------------
  intrea = 0
  do i = 1, CARRAY_MXCOLM
     if(cwork(2)(i:i) == '.') then
        intrea = 1
     end if
  end do

  if(intrea <= 0) then
     error_message = 'Error: in reading: ' // trim(cvalue)
     call util_error(ERROR%STOP_STD, error_message)
  end if

  !---- get value --------------------------------------------------------
  read (cwork(2), *) rvalue
  write (lunout, '(3a, g10.3)') '---reading real parameter: ', trim(cvalue), ' = ', rvalue

end subroutine ukoto_rvalue2


subroutine ukoto_ivalue2(lunout, cwork, ivalue, cvalue)

  use const_maxsize
  use const_index
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lunout
  integer, intent(out) :: ivalue
  character(CARRAY_MXCOLM), intent(in) :: cwork(2)
  character(CARRAY_MXCOLM), intent(in) :: cvalue

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i, len
  integer :: intrea
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
!  len = len_trim(cwork(1))

!  if(cwork(1)(1:len) /= cvalue(1:len)) then
  if(cwork(1) /= cvalue) then
     return
  end if

  !---- store characters into 'ctmp01' -----------------------------------
  intrea = 0
  do i = 1, CARRAY_MXCOLM
     if(cwork(2)(i:i) == '.') then
        intrea = 1
     end if
  end do

  if(intrea > 0) then
     error_message = 'Error: in reading ' // trim(cvalue)
     call util_error(ERROR%STOP_STD, error_message)
  end if

  !---- get value --------------------------------------------------------
  read (cwork(2), *) ivalue
  len = len_trim(cwork(2))
  if(len <= 1) then
     write (lunout, '(3a, i1)') '---reading integer parameter: ', trim(cvalue), ' = ', ivalue
  else if(len <= 4) then
     write (lunout, '(3a, i4)') '---reading integer parameter: ', trim(cvalue), ' = ', ivalue
  else
     write (lunout, '(3a, i10)') '---reading integer parameter: ', trim(cvalue), ' = ', ivalue
  end if
end subroutine ukoto_ivalue2


subroutine ukoto_lvalue2(lunout, cwork, lvalue, cvalue)

  use const_maxsize
  use const_index
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lunout
  integer(L_INT), intent(out) :: lvalue
  character(CARRAY_MXCOLM), intent(in) :: cwork(2)
  character(CARRAY_MXCOLM), intent(in) :: cvalue

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i  !, len
  integer :: intrea
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
!  len = len_trim(cwork(1))

!  if(cwork(1)(1:len) /= cvalue(1:len)) then
  if(cwork(1) /= cvalue) then
     return
  end if

  !---- store characters into 'ctmp01' -----------------------------------
  intrea = 0
  do i = 1, CARRAY_MXCOLM
     if(cwork(2)(i:i) == '.') then
        intrea = 1
     end if
  end do

  if(intrea > 0) then
     error_message = 'Error: in reading ' // trim(cvalue)
     call util_error(ERROR%STOP_STD, error_message)
  end if

  !---- get value --------------------------------------------------------
  read (cwork(2), *) lvalue
  write (lunout, '(3a, i12)') '---reading long parameter: ', trim(cvalue), ' = ', lvalue

end subroutine ukoto_lvalue2

